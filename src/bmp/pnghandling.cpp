#ifdef HAVE_PNG

#include "pnghandling.h"
#include <assert.h>

PngHandling::PngHandling() : width(1),height(1), color_type(PNG_COLOR_TYPE_RGB_ALPHA), bit_depth(8)
{
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_ptr = png_create_info_struct(png_ptr);

    if (!info_ptr)
        cerr <<  "Error creating png_create_info_struct failed" << endl;
}

PngHandling::~PngHandling()
{
    /* cleanup heap allocation */
    for (int y=0; y<height; y++){
        free(row_pointers[y]);
    }
    free(row_pointers);

    png_destroy_write_struct(&png_ptr,&info_ptr);

}

bool PngHandling::isInitialized() const
{
    if (!png_ptr)   {
        cerr <<  "Error: creating png_ptr failed" << endl;
        return false;
    }
    if (!info_ptr){
        cerr <<  "Error: creating png_create_info_struct failed" << endl;
        return false;
    }
    return true;
}

void PngHandling::setImageDimensions(uint w, uint h)
{
    width=w;
    height=h;
}

void PngHandling::prepareImage()
{
    //allocate png byte** 2d array containing image data

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) malloc( sizeof(png_byte)*4*width);

    //putBackground(veryBlack);
}

void PngHandling::putPixel(uint x, uint y, const Color& c)
{
    assert(isInitialized());

    if (!(x<width and y<height))
        return;

    png_byte* row = row_pointers[y];
    png_byte* ptr = &(row[x*4]);

    unsigned int p0,p1,p2;
    c.getIntegerRGB(p0,p1,p2);
    ptr[0]=p0; ptr[1]=p1; ptr[2]=p2; ptr[3]=255;
}

void PngHandling::writeToFile(char* file_name)
{
    assert(isInitialized());
    /* create file */
    FILE *fp = fopen(file_name, "wb");
    if (!fp)
            cerr << "Error opening file: " << string(file_name) << endl;

    //FIXME After Ubuntu update setjmp stopped returned errors so temporarily removed

    //set up error handling with setjmp
    /*if (setjmp(png_jmpbuf(png_ptr))){
        png_destroy_write_struct(&png_ptr,&info_ptr);
        fclose(fp);
        cerr << "Error during init_io" << endl;;
    }*/

    png_init_io(png_ptr, fp);

    //Setting the contents of info for output
    /*if (setjmp(png_jmpbuf(png_ptr)))
        cerr << "Error during writing header" << endl;;
    */
    png_set_IHDR(png_ptr, info_ptr, width, height,
            bit_depth, color_type, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);


    /* write bytes */
    /*if (setjmp(png_jmpbuf(png_ptr)))
        cerr << "Error during writing bytes" << endl;
    */
    png_write_image(png_ptr, row_pointers);

    /*if (setjmp(png_jmpbuf(png_ptr)))
        cerr << "Error during end of write" << endl;
    */
    png_write_end(png_ptr, NULL);
    fclose(fp);
}

void PngHandling::putBackground(const Color& c)
{
    for (int y=0; y<height; y++) {
        png_byte* row = row_pointers[y];
        for (int x=0; x<width; x++) {
            png_byte* ptr = &(row[x*4]);
            uint p0,p1,p2;
            c.getIntegerRGB(p0,p1,p2);
            ptr[0]=p0; ptr[1]=p1; ptr[2]=p2; ptr[3]=255;
        }
    }
}

#endif
