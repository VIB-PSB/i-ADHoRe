#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include "DataSet.h"
#include "util.h"
#include "parallel.h"
#include "Settings.h"
#include <pthread.h>

extern "C" void* startThread(void *args)
{
    ThreadArgs *threadArgs = reinterpret_cast<ThreadArgs*>(args);
    DataSet *dataset = threadArgs->dataset;
    int threadID = threadArgs->threadID;

    vector<Multiplicon *> threadMultiplicons;
    vector<SynthenicCloud*> threadClouds;

    while (true) {
        // perform all the work there is to perform
        int firstItem, nItems; bool threadDidWork = false;
        while (dataset->runGetSomeWork(firstItem, nItems,
                                       threadID, !threadDidWork)) {
            dataset->runWorkFunction(firstItem, nItems, threadID, threadMultiplicons, threadClouds);
            threadDidWork = true;
        }

        // if the thread did no work, unnecessary to lock mutexes etc...
        if (threadDidWork) {
            // protect again simulatenous writing
            pthread_mutex_lock (&dataset->mpcMutex);
            dataset->localMultiplicons.insert(dataset->localMultiplicons.end(),
                                              threadMultiplicons.begin(),
                                              threadMultiplicons.end());
            dataset->localClouds.insert(dataset->localClouds.end(),
                                        threadClouds.begin(),
                                        threadClouds.end());
            pthread_mutex_unlock (&dataset->mpcMutex);

            // clear the found multiplicons for the next run
            threadMultiplicons.clear();
            threadClouds.clear();

            pthread_mutex_lock (&dataset->queueMutex);
            dataset->workInProgress--;
            // signal master thread that threads are finished
            if (dataset->workInProgress == 0)
                pthread_cond_signal (&dataset->masterCond);
            pthread_mutex_unlock (&dataset->queueMutex);
        }

        pthread_mutex_lock (&dataset->workerMutex);
        if (dataset->destroyTP) {
            pthread_mutex_unlock (&dataset->workerMutex);
            pthread_exit(NULL);
        }

        // go to sleep
        // a spurious wakeup is not that bad here, the thread will wake up,
        // find out that there is no work, and go back to sleep
        pthread_cond_wait( &dataset->workerCond, &dataset->workerMutex );
        pthread_mutex_unlock (&dataset->workerMutex);
    }

    pthread_exit(NULL);
}

void DataSet::wakeThreads()
{
    // wake up the slaves
    pthread_mutex_lock(&workerMutex);
    pthread_cond_broadcast(&workerCond);
    pthread_mutex_unlock(&workerMutex);
}

void DataSet::finishWorkerThreads()
{
    // make sure your worker threads are finished before continuing
    pthread_mutex_lock(&queueMutex);
    while (workInProgress != 0)
        pthread_cond_wait(&masterCond, &queueMutex);
    pthread_mutex_unlock(&queueMutex);
}

void DataSet::createThreadPool()
{
    pthread_mutex_init(&mpcMutex, NULL);
    pthread_mutex_init(&workerMutex, NULL);
    pthread_mutex_init(&queueMutex, NULL);
    pthread_cond_init(&workerCond, NULL);
    pthread_cond_init(&masterCond, NULL);

    // spawn extra threads, if necessary
    nThreads = settings.getNumThreads() - 1;

    if (nThreads > 0) {
        threads = new pthread_t[nThreads];
        threadArgs = new ThreadArgs[nThreads];

        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

        // create the worker threads without any workload
        destroyTP = false;
        currentIndex = 1; finalIndex = 0;
        for (int threadID = 0; threadID < nThreads; threadID++) {
            threadArgs[threadID].dataset = this;
            threadArgs[threadID].threadID = threadID + 1; // 0 == master thread
            if (pthread_create(&threads[threadID], &attr,
                               *startThread, &threadArgs[threadID]) != 0) {
                cerr << "Failed to create thread " << threadID << "/" << nThreads << endl;
                cerr << "Try using fewer threads..." << endl;
                exit(EXIT_FAILURE);
            }
        }

        pthread_attr_destroy(&attr);
    }
}

void DataSet::destroyThreadPool()
{
    if (nThreads > 0) {
        // signal the worker threads that it is time to terminate
        pthread_mutex_lock (&workerMutex);
        destroyTP = true;
        pthread_cond_broadcast (&workerCond);
        pthread_mutex_unlock (&workerMutex);

        // wait for the threads to terminate
        void *status;
        for (int threadID = 0; threadID < nThreads; threadID++)
            pthread_join(threads[threadID], &status);

        // destroy all resources
        delete [] threads;
        delete [] threadArgs;
        threads = NULL;
        threadArgs = NULL;
    }

    pthread_mutex_destroy(&mpcMutex);
    pthread_mutex_destroy(&workerMutex);
    pthread_mutex_destroy(&queueMutex);
    pthread_cond_destroy(&workerCond);
    pthread_cond_destroy(&masterCond);
}
