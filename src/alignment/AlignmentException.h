#ifndef __ALIGNMENTEXCEPTION_H
#define __ALIGNMENTEXCEPTION_H


class AlignmentException : public runtime_error {
    public: AlignmentException(const string& msg = "") : runtime_error(msg) {}
};


#endif
