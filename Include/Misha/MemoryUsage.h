/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef MEMORY_USAGE_INCLUDED
#define MEMORY_USAGE_INCLUDED

#ifdef WIN32

#include <Windows.h>
#include <psapi.h>

class MemoryInfo
{
public:
    size_t TotalPhysicalMemory;
    size_t FreePhysicalMemory;
    size_t TotalSwapSpace;
    size_t FreeSwapSpace;
    size_t TotalVirtualAddressSpace;
    size_t FreeVirtualAddressSpace;
    size_t PageSize;

    void set(void)
    {
        MEMORYSTATUSEX Mem;
        SYSTEM_INFO Info;
        ZeroMemory( &Mem, sizeof(Mem));
        ZeroMemory( &Info, sizeof(Info)); 
        Mem.dwLength = sizeof(Mem);
        ::GlobalMemoryStatusEx( &Mem );
        ::GetSystemInfo( &Info );

        TotalPhysicalMemory = (size_t)Mem.ullTotalPhys;
        FreePhysicalMemory = (size_t)Mem.ullAvailPhys;
        TotalSwapSpace = (size_t)Mem.ullTotalPageFile;
        FreeSwapSpace = (size_t)Mem.ullAvailPageFile;
        TotalVirtualAddressSpace = (size_t)Mem.ullTotalVirtual;
        FreeVirtualAddressSpace = (size_t)Mem.ullAvailVirtual;
        PageSize = (size_t)Info.dwPageSize;
    }
    size_t usage(void) const {return TotalVirtualAddressSpace-FreeVirtualAddressSpace;}

    static size_t Usage(void){
        MEMORY_BASIC_INFORMATION mbi; 
        size_t      dwMemUsed = 0; 
        PVOID      pvAddress = 0; 


        memset(&mbi, 0, sizeof(MEMORY_BASIC_INFORMATION)); 
        while(VirtualQuery(pvAddress, &mbi, sizeof(MEMORY_BASIC_INFORMATION)) == sizeof(MEMORY_BASIC_INFORMATION)){ 
            if(mbi.State == MEM_COMMIT && mbi.Type == MEM_PRIVATE){dwMemUsed += mbi.RegionSize;}
            pvAddress = ((BYTE*)mbi.BaseAddress) + mbi.RegionSize; 
        } 
        return dwMemUsed; 
    }

    ////////////////////////////////////////////////////////////////////////////
    /*! Determines the amount of used memory in megabytes.
    //  @return amount of resident memory (megabytes)
    *///////////////////////////////////////////////////////////////////////////
    static double UsageMB(void)
    {
        double u = ((double) Usage()) / (1 << 20);
        if( u>peakMB ) peakMB = u;
        return u;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Gets the maximum megabytes of memory used during the program's run
    ////////////////////////////////////////////////////////////////////////////
    static size_t PeakUsage( void )
    {
        HANDLE hProcess = GetCurrentProcess();
        PROCESS_MEMORY_COUNTERS pmc;

        if ( GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc)) ) return  pmc.PeakWorkingSetSize;
        else fprintf( stderr , "Failed to get peak working size\n" ) , exit( 0 );
    }
    ////////////////////////////////////////////////////////////////////////////
    // Gets the maximum megabytes of memory used during the program's run
    ////////////////////////////////////////////////////////////////////////////
    static double PeakUsageMB( void )
    {
		return ((double)PeakUsage()) / (1<<20);
    }

    static double peakMB;
};
double MemoryInfo::peakMB = 0;
#else   // !WIN32

#ifdef __APPLE__    // !__APPLE__
#include <sys/resource.h>
#include <sys/sysctl.h>
#include <mach/mach.h>
#include <cstdio>
#include <stdint.h>
#include <unistd.h>

////////////////////////////////////////////////////////////////////////////////
/*! Class for retrieving memory utilization on Mac OS X.
//  (Or probably any system using the MACH kernel.)
*///////////////////////////////////////////////////////////////////////////////
class MemoryInfo
{
    public:
        static double peakMB;
        /** The amount of physical memory (bytes) */
        uint64_t mem_physical;

        ////////////////////////////////////////////////////////////////////////
        /*! Constructor queries various memory properties of the system
        *///////////////////////////////////////////////////////////////////////
        MemoryInfo()
        {
            int mib[2];
            mib[0] = CTL_HW;
            mib[1] = HW_MEMSIZE;

            size_t returnSize = sizeof(mem_physical);
            if (sysctl(mib, 2, &mem_physical, &returnSize, NULL, 0) == -1)
                perror("Error in sysctl call");
        }

        ////////////////////////////////////////////////////////////////////////
        /*! Queries the kernel for the amount of resident memory in bytes.
        //  @return amount of resident memory (bytes)
        *///////////////////////////////////////////////////////////////////////
        static size_t Usage(void)
        {
            task_t targetTask = mach_task_self();
            struct task_basic_info ti;
            mach_msg_type_number_t count = TASK_BASIC_INFO_64_COUNT;

            kern_return_t kr = task_info(targetTask, TASK_BASIC_INFO_64,
                                         (task_info_t) &ti, &count);
            if (kr != KERN_SUCCESS) {
                printf("Kernel returned error during memory usage query");
                return -1;
            }

            // On Mac oS X, the resident_size is in bytes, not pages!
            // (This is different from the GNU kernel)
            return ti.resident_size;
        }

        ////////////////////////////////////////////////////////////////////////
        /*! Queries the kernel for the amount of resident memory in megabytes.
        //  @return amount of resident memory (megabytes)
        *///////////////////////////////////////////////////////////////////////
        static double UsageMB(void)
        {
            double u = ((double) Usage()) / (1 << 20);
            if( u>peakMB ) peakMB = u;
            return u;
        }

        ////////////////////////////////////////////////////////////////////////
        // Gets the maximum megabytes of memory used during the program's run
        ////////////////////////////////////////////////////////////////////////
        static double PeakUsageMB(void)
        {
            struct rusage usageInfo;
            getrusage(RUSAGE_SELF, &usageInfo);
            double peakKB = usageInfo.ru_maxrss;
            return peakKB / (1 << 20);
        }
};
double MemoryInfo::peakMB = 0;

#else   // !__APPLE__

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

class MemoryInfo
{
    public:
        static double peakMB;
        ////////////////////////////////////////////////////////////////////////
        /*! Gets the amount of resident memory (in bytes)
        //  @return     resident memory (in bytes)
        *///////////////////////////////////////////////////////////////////////
        static size_t Usage(void)
        {
            unsigned long vm = 0;
            FILE* f = fopen("/proc/self/stat","rb");

            long rss;
            if (f) {
                int n = fscanf(f, "%*d %*s %*c %*d %*d %*d %*d %*d %*lu %*lu %*lu %*lu %*lu %*lu %*lu %*ld %*ld %*ld %*ld %*ld %*ld %*lu %*lu %ld", &rss);
                fclose(f);
            }
            return rss * getpagesize();
        }

        ////////////////////////////////////////////////////////////////////////////
        /*! Determines the amount of used memory in megabytes.
        //  @return amount of resident memory (megabytes)
        *///////////////////////////////////////////////////////////////////////////
        static double UsageMB(void)
        {
            double u = ((double) Usage()) / (1 << 20);
            if (u>peakMB)   peakMB = u;
            return u;
        }

        ////////////////////////////////////////////////////////////////////////
        // Gets the maximum megabytes of memory used during the program's run
        ////////////////////////////////////////////////////////////////////////
        static double PeakUsageMB(void)
        {
            struct rusage usageInfo;
            getrusage(RUSAGE_SELF, &usageInfo);
            double peakKB = usageInfo.ru_maxrss;
            return peakKB / (1 << 10);
        }
};
double MemoryInfo::peakMB = 0;
#endif // __APPLE__
#endif // WIN32

#endif // MEMORY_USAGE_INCLUDE
