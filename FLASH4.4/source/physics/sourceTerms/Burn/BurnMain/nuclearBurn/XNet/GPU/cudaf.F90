MODULE cudaf
!===============================================================================
! INTERFACE to CUDA C routines
!===============================================================================
    USE, INTRINSIC :: ISO_C_BINDING

    ENUM, BIND(C) !:: cudaMemcpyKind
    ENUMERATOR :: cudaMemcpyHostToHost=0
    ENUMERATOR :: cudaMemcpyHostToDevice=1
    ENUMERATOR :: cudaMemcpyDeviceToHost=2
    ENUMERATOR :: cudaMemcpyDeviceToDevice=3
    ENUMERATOR :: cudaMemcpyDefault=4
    END ENUM !cudaMemcpyKind

    ENUM, BIND(C) !:: cudaHostFlags
    ENUMERATOR :: cudaHostAllocDefault=INT(z'00')
    ENUMERATOR :: cudaHostAllocPortable=INT(z'01')
    ENUMERATOR :: cudaHostAllocMapped=INT(z'02')
    ENUMERATOR :: cudaHostAllocWriteCombined=INT(z'04')
    END ENUM !cudaHostFlags

    ENUM, BIND(C) !:: cudaDeviceFlags
    ENUMERATOR :: cudaDeviceScheduleAuto=INT(z'00')
    ENUMERATOR :: cudaDeviceScheduleSpin=INT(z'01')
    ENUMERATOR :: cudaDeviceScheduleYield=INT(z'02')
    ENUMERATOR :: cudaDeviceScheduleBlockingSync=INT(z'04')
    ENUMERATOR :: cudaDeviceBlockingSync=INT(z'04')
    ENUMERATOR :: cudaDeviceScheduleMask=INT(z'07')
    ENUMERATOR :: cudaDeviceMapHost=INT(z'08')
    ENUMERATOR :: cudaDeviceLmemResizeToMax=INT(z'10')
    ENUMERATOR :: cudaDeviceMask=INT(z'1f')
    END ENUM !cudaDeviceFlags

    ENUM, BIND(C) !:: cudaStreamFlags
    ENUMERATOR :: cudaStreamDefault=INT(z'00')
    ENUMERATOR :: cudaStreamNonBlocking=INT(z'01')
    END ENUM

    ENUM, BIND(C) !:: cudaEventFlags
    ENUMERATOR :: cudaEventDefault=INT(z'00')
    ENUMERATOR :: cudaEventBlockingSync=INT(z'01')
    ENUMERATOR :: cudaEventDisableTiming=INT(z'02')
    ENUMERATOR :: cudaEventInterprocess=INT(z'04')
    END ENUM

    ENUM, BIND(C) !:: cudaSharedMemConfig
    ENUMERATOR :: cudaSharedMemBankSizeDefault   = 0
    ENUMERATOR :: cudaSharedMemBankSizeFourByte  = 1
    ENUMERATOR :: cudaSharedMemBankSizeEightByte = 2
    END ENUM

    ENUM, BIND(C) !:: cudaFuncCache
    ENUMERATOR :: cudaFuncCachePreferNone   = 0
    ENUMERATOR :: cudaFuncCachePreferShared = 1
    ENUMERATOR :: cudaFuncCachePreferL1     = 2
    ENUMERATOR :: cudaFuncCachePreferEqual  = 3
    END ENUM

    ENUM, BIND(C) !:: cudaError
    ENUMERATOR :: cudaSuccess                           =      0
    ENUMERATOR :: cudaErrorMissingConfiguration         =      1
    ENUMERATOR :: cudaErrorMemoryAllocation             =      2
    ENUMERATOR :: cudaErrorInitializationError          =      3
    ENUMERATOR :: cudaErrorLaunchFailure                =      4
    ENUMERATOR :: cudaErrorPriorLaunchFailure           =      5
    ENUMERATOR :: cudaErrorLaunchTimeout                =      6
    ENUMERATOR :: cudaErrorLaunchOutOfResources         =      7
    ENUMERATOR :: cudaErrorInvalidDeviceFunction        =      8
    ENUMERATOR :: cudaErrorInvalidConfiguration         =      9
    ENUMERATOR :: cudaErrorInvalidDevice                =     10
    ENUMERATOR :: cudaErrorInvalidValue                 =     11
    ENUMERATOR :: cudaErrorInvalidPitchValue            =     12
    ENUMERATOR :: cudaErrorInvalidSymbol                =     13
    ENUMERATOR :: cudaErrorMapBufferObjectFailed        =     14
    ENUMERATOR :: cudaErrorUnmapBufferObjectFailed      =     15
    ENUMERATOR :: cudaErrorInvalidHostPointer           =     16
    ENUMERATOR :: cudaErrorInvalidDevicePointer         =     17
    ENUMERATOR :: cudaErrorInvalidTexture               =     18
    ENUMERATOR :: cudaErrorInvalidTextureBinding        =     19
    ENUMERATOR :: cudaErrorInvalidChannelDescriptor     =     20
    ENUMERATOR :: cudaErrorInvalidMemcpyDirection       =     21
    ENUMERATOR :: cudaErrorAddressOfConstant            =     22
    ENUMERATOR :: cudaErrorTextureFetchFailed           =     23
    ENUMERATOR :: cudaErrorTextureNotBound              =     24
    ENUMERATOR :: cudaErrorSynchronizationError         =     25
    ENUMERATOR :: cudaErrorInvalidFilterSetting         =     26
    ENUMERATOR :: cudaErrorInvalidNormSetting           =     27
    ENUMERATOR :: cudaErrorMixedDeviceExecution         =     28
    ENUMERATOR :: cudaErrorCudartUnloading              =     29
    ENUMERATOR :: cudaErrorUnknown                      =     30
    ENUMERATOR :: cudaErrorNotYetImplemented            =     31
    ENUMERATOR :: cudaErrorMemoryValueTooLarge          =     32
    ENUMERATOR :: cudaErrorInvalidResourceHandle        =     33
    ENUMERATOR :: cudaErrorNotReady                     =     34
    ENUMERATOR :: cudaErrorInsufficientDriver           =     35
    ENUMERATOR :: cudaErrorSetOnActiveProcess           =     36
    ENUMERATOR :: cudaErrorInvalidSurface               =     37
    ENUMERATOR :: cudaErrorNoDevice                     =     38
    ENUMERATOR :: cudaErrorECCUncorrectable             =     39
    ENUMERATOR :: cudaErrorSharedObjectSymbolNotFound   =     40
    ENUMERATOR :: cudaErrorSharedObjectInitFailed       =     41
    ENUMERATOR :: cudaErrorUnsupportedLimit             =     42
    ENUMERATOR :: cudaErrorDuplicateVariableName        =     43
    ENUMERATOR :: cudaErrorDuplicateTextureName         =     44
    ENUMERATOR :: cudaErrorDuplicateSurfaceName         =     45
    ENUMERATOR :: cudaErrorDevicesUnavailable           =     46
    ENUMERATOR :: cudaErrorInvalidKernelImage           =     47
    ENUMERATOR :: cudaErrorNoKernelImageForDevice       =     48
    ENUMERATOR :: cudaErrorIncompatibleDriverContext    =     49
    ENUMERATOR :: cudaErrorPeerAccessAlreadyEnabled     =     50
    ENUMERATOR :: cudaErrorPeerAccessNotEnabled         =     51
    ENUMERATOR :: cudaErrorDeviceAlreadyInUse           =     54
    ENUMERATOR :: cudaErrorProfilerDisabled             =     55
    ENUMERATOR :: cudaErrorProfilerNotInitialized       =     56
    ENUMERATOR :: cudaErrorProfilerAlreadyStarted       =     57
    ENUMERATOR :: cudaErrorProfilerAlreadyStopped       =     58
    ENUMERATOR :: cudaErrorAssert                       =     59
    ENUMERATOR :: cudaErrorTooManyPeers                 =     60
    ENUMERATOR :: cudaErrorHostMemoryAlreadyRegistered  =     61
    ENUMERATOR :: cudaErrorHostMemoryNotRegistered      =     62
    ENUMERATOR :: cudaErrorOperatingSystem              =     63
    ENUMERATOR :: cudaErrorPeerAccessUnsupported        =     64
    ENUMERATOR :: cudaErrorLaunchMaxDepthExceeded       =     65
    ENUMERATOR :: cudaErrorLaunchFileScopedTex          =     66
    ENUMERATOR :: cudaErrorLaunchFileScopedSurf         =     67
    ENUMERATOR :: cudaErrorSyncDepthExceeded            =     68
    ENUMERATOR :: cudaErrorLaunchPendingCountExceeded   =     69
    ENUMERATOR :: cudaErrorNotPermitted                 =     70
    ENUMERATOR :: cudaErrorNotSupported                 =     71
    ENUMERATOR :: cudaErrorStartupFailure               =   INT(z'7f')
    ENUMERATOR :: cudaErrorApiFailureBase               =  10000
    END ENUM !cudaError

    TYPE, bind(C) :: cudaDeviceProp
        Character(C_CHAR) :: name(256)
        INTEGER(C_SIZE_T) :: totalGlobalMem
        INTEGER(C_SIZE_T) :: sharedMemPerBlock
        INTEGER(C_INT) :: regsPerBlock
        INTEGER(C_INT) :: warpSize
        INTEGER(C_SIZE_T) :: memPitch
        INTEGER(C_INT) :: maxThreadsPerBlock
        INTEGER(C_INT) :: maxThreadsDim(3)
        INTEGER(C_INT) :: maxGridSize(3)
        INTEGER(C_INT) :: clockRate
        INTEGER(C_SIZE_T) :: totalConstMem
        INTEGER(C_INT) :: major
        INTEGER(C_INT) :: minor
        INTEGER(C_SIZE_T) :: textureAlignment
        INTEGER(C_SIZE_T) :: texturePitchAlignment
        INTEGER(C_INT) :: deviceOverlap
        INTEGER(C_INT) :: multiProcessorCount
        INTEGER(C_INT) :: kernelExecTimeoutEnabled
        INTEGER(C_INT) :: integrated
        INTEGER(C_INT) :: canMapHostMemory
        INTEGER(C_INT) :: computeMode
        INTEGER(C_INT) :: maxTexture1D
        INTEGER(C_INT) :: maxTexture1DMipmap
        INTEGER(C_INT) :: maxTexture1DLinear
        INTEGER(C_INT) :: maxTexture2D(2)
        INTEGER(C_INT) :: maxTexture2DMipmap(2)
        INTEGER(C_INT) :: maxTexture2DLinear(2)
        INTEGER(C_INT) :: maxTexture2DGather(2)
        INTEGER(C_INT) :: maxTexture3D(3)
        INTEGER(C_INT) :: maxTexture3DAlt(3)
        INTEGER(C_INT) :: maxTextureCubemap
        INTEGER(C_INT) :: maxTexture1DLayered(2)
        INTEGER(C_INT) :: maxTexture2DLayered(3)
        INTEGER(C_INT) :: maxTextureCubemapLayered(2)
        INTEGER(C_INT) :: maxSurface1D
        INTEGER(C_INT) :: maxSurface2D(2)
        INTEGER(C_INT) :: maxSurface3D(3)
        INTEGER(C_INT) :: maxSurface1DLayered(2)
        INTEGER(C_INT) :: maxSurface2DLayered(3)
        INTEGER(C_INT) :: maxSurfaceCubemap
        INTEGER(C_INT) :: maxSurfaceCubemapLayered(2)
        INTEGER(C_SIZE_T) :: surfaceAlignment
        INTEGER(C_INT) :: concurrentKernels
        INTEGER(C_INT) :: ECCEnabled
        INTEGER(C_INT) :: pciBusID
        INTEGER(C_INT) :: pciDeviceID
        INTEGER(C_INT) :: pciDomainID
        INTEGER(C_INT) :: tccDriver
        INTEGER(C_INT) :: asyncEngineCount
        INTEGER(C_INT) :: unifiedAddressing
        INTEGER(C_INT) :: memoryClockRate
        INTEGER(C_INT) :: memoryBusWidth
        INTEGER(C_INT) :: l2CacheSize
        INTEGER(C_INT) :: maxThreadsPerMultiProcessor
        INTEGER(C_INT) :: streamPrioritiesSupported
        INTEGER(C_INT) :: globalL1CacheSupported
        INTEGER(C_INT) :: localL1CacheSupported
        INTEGER(C_SIZE_T) :: sharedMemPerMultiprocessor
        INTEGER(C_INT) :: regsPerMultiprocessor
        INTEGER(C_INT) :: managedMemSupported
        INTEGER(C_INT) :: isMultiGpuBoard
        INTEGER(C_INT) :: multiGpuBoardGroupID
        INTEGER(C_INT) :: singleToDoublePrecisionPerfRatio
        INTEGER(C_INT) :: pageableMemoryAccess
        INTEGER(C_INT) :: concurrentManagedAccess
    END TYPE cudaDeviceProp

    INTERFACE

        FUNCTION cudaHostAlloc(cPtr, size, flags) &
        &   BIND(C, NAME="cudaHostAlloc")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaHostAlloc
            TYPE(C_PTR) :: cPtr
            INTEGER(C_SIZE_T), VALUE :: size
            INTEGER(C_INT), VALUE :: flags
        END FUNCTION cudaHostAlloc

        FUNCTION cudaMallocHost(cPtr, size) &
        &   BIND(C, NAME="cudaMallocHost")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaMallocHost
            TYPE(C_PTR) :: cPtr
            INTEGER(C_SIZE_T), VALUE :: size
        END FUNCTION cudaMallocHost

        FUNCTION cudaFreeHost(cPtr) &
        &   BIND(C, NAME="cudaFreeHost")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaFreeHost
            TYPE(C_PTR), VALUE :: cPtr
        END FUNCTION cudaFreeHost

        FUNCTION cudaMalloc(dPtr, size) &
        &   BIND(C, NAME="cudaMalloc")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaMalloc
            TYPE(C_PTR) :: dPtr
            INTEGER(C_SIZE_T), VALUE :: size
        END FUNCTION cudaMalloc

        FUNCTION cudaFree(dPtr) &
        &   BIND(C, NAME="cudaFree")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaFree
            TYPE(C_PTR), VALUE :: dPtr
        END FUNCTION cudaFree

        FUNCTION cudaMemcpy(dst, src, memSize, cpyKind) &
        &   BIND(C, NAME="cudaMemcpy")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaMemcpy
            TYPE(C_PTR), VALUE :: dst
            TYPE(C_PTR), VALUE :: src
            INTEGER(C_SIZE_T), VALUE :: memSize
            INTEGER(kind(cudaMemcpyHostToHost)), VALUE :: cpyKind
        END FUNCTION cudaMemcpy

        FUNCTION cudaMemcpyAsync(dst, src, memSize, cpyKind, stream) &
        &   BIND(C, NAME="cudaMemcpyAsync")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaMemcpyAsync
            TYPE(C_PTR), VALUE :: dst
            TYPE(C_PTR), VALUE :: src
            INTEGER(C_SIZE_T), VALUE :: memSize
            INTEGER(kind(cudaMemcpyHostToHost)), VALUE :: cpyKind
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cudaMemcpyAsync

        FUNCTION cudaStreamCreate(stream) &
        &   BIND(C, NAME="cudaStreamCreate")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaStreamCreate
            TYPE(C_PTR) :: stream
        END FUNCTION cudaStreamCreate

        FUNCTION cudaStreamCreateWithFlags(stream, flags) &
        &   BIND(C, NAME="cudaStreamCreateWithFlags")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaStreamCreateWithFlags
            TYPE(C_PTR) :: stream
            INTEGER(C_INT), VALUE :: flags
        END FUNCTION cudaStreamCreateWithFlags

        FUNCTION cudaStreamDestroy(stream) &
        &   BIND(C, NAME="cudaStreamDestroy")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaStreamDestroy
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cudaStreamDestroy

        FUNCTION cudaStreamSynchronize(stream) &
        &   BIND(C, NAME="cudaStreamSynchronize")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaStreamSynchronize
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cudaStreamSynchronize

        FUNCTION cudaStreamWaitEvent(stream, event) &
        &   BIND(C, NAME="cudaStreamWaitEvent")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaStreamWaitEvent
            TYPE(C_PTR), VALUE :: stream
            TYPE(C_PTR), VALUE :: event
        END FUNCTION cudaStreamWaitEvent

        FUNCTION cudaGetDeviceCount(deviceCount) &
        &   BIND(C, NAME="cudaGetDeviceCount")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaGetDeviceCount
            INTEGER(C_INT) :: deviceCount
        END FUNCTION cudaGetDeviceCount

        FUNCTION cudaGetDeviceProperties(prop, device) &
        &   BIND(C, NAME="cudaGetDeviceProperties")
            USE, INTRINSIC :: ISO_C_BINDING
            IMPORT cudaDeviceProp
            INTEGER(KIND(cudaSuccess)) :: cudaGetDeviceProperties
            TYPE(cudaDeviceProp) :: prop
            INTEGER(C_INT), VALUE :: device
        END FUNCTION cudaGetDeviceProperties

        FUNCTION cudaSetDevice(device) &
        &   BIND(C, NAME="cudaSetDevice")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaSetDevice
            INTEGER(C_INT), VALUE :: device
        END FUNCTION cudaSetDevice

        FUNCTION cudaSetDeviceFlags(flags) &
        &   BIND(C, NAME="cudaSetDeviceFlags")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaSetDeviceFlags
            INTEGER(C_INT), VALUE :: flags
        END FUNCTION cudaSetDeviceFlags

        FUNCTION cudaHostGetDevicePointer(dPtr, cPtr, flags) &
        &   BIND(C, NAME="cudaHostGetDevicePointer")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaHostGetDevicePointer
            TYPE(C_PTR) :: dPtr
            TYPE(C_PTR), VALUE :: cPtr
            INTEGER(C_INT), VALUE :: flags
        END FUNCTION cudaHostGetDevicePointer

        FUNCTION cudaDeviceSynchronize() &
        &   BIND(C, NAME="cudaDeviceSynchronize")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaDeviceSynchronize
        END FUNCTION cudaDeviceSynchronize

        FUNCTION cudaDeviceReset() &
        &   BIND(C, NAME="cudaDeviceReset")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaDeviceReset
        END FUNCTION cudaDeviceReset

        FUNCTION cudaEventCreate( event ) &
        &   BIND(C, NAME="cudaEventCreate")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaEventCreate
            TYPE(C_PTR) :: event
        END FUNCTION cudaEventCreate

        FUNCTION cudaEventCreateWithFlags( event, flags ) &
        &   BIND(C, NAME="cudaEventCreateWithFlags")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaEventCreateWithFlags
            TYPE(C_PTR) :: event
            INTEGER(C_INT), VALUE :: flags
        END FUNCTION cudaEventCreateWithFlags

        FUNCTION cudaEventDestroy( event ) &
        &   BIND(C, NAME="cudaEventDestroy")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaEventDestroy
            TYPE(C_PTR), VALUE :: event
        END FUNCTION cudaEventDestroy

        FUNCTION cudaEventElapsedTime( eventTime, eventStart, eventStop ) &
        &   BIND(C, NAME="cudaEventElapsedTime")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaEventElapsedTime
            REAL(C_FLOAT) :: eventTime
            TYPE(C_PTR), VALUE :: eventStart
            TYPE(C_PTR), VALUE :: eventStop
        END FUNCTION cudaEventElapsedTime

        FUNCTION cudaEventQuery( event ) &
        &   BIND(C, NAME="cudaEventQuery")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaEventQuery
            TYPE(C_PTR), VALUE :: event
        END FUNCTION cudaEventQuery

        FUNCTION cudaEventRecord( event, stream ) &
        &   BIND(C, NAME="cudaEventRecord")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaEventRecord
            TYPE(C_PTR), VALUE :: event
            TYPE(C_PTR), VALUE :: stream
        END FUNCTION cudaEventRecord

        FUNCTION cudaEventSynchronize( event ) &
        &   BIND(C, NAME="cudaEventSynchronize")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaEventSynchronize
            TYPE(C_PTR), VALUE :: event
        END FUNCTION cudaEventSynchronize

        FUNCTION cudaDeviceSetSharedMemConfig( config ) &
        &   BIND(C, NAME="cudaDeviceSetSharedMemConfig")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaDeviceSetSharedMemConfig
            INTEGER(C_INT), VALUE :: config
        END FUNCTION cudaDeviceSetSharedMemConfig

        FUNCTION cudaDeviceSetCacheConfig( cacheConfig ) &
        &   BIND(C, NAME="cudaDeviceSetCacheConfig")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND(cudaSuccess)) :: cudaDeviceSetCacheConfig
            INTEGER(C_INT), VALUE :: cacheConfig
        END FUNCTION cudaDeviceSetCacheConfig

    END INTERFACE

END MODULE cudaf
