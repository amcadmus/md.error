/**
 * @file   NeighborList_interface.cu
 * @author Han Wang <han_wang@math.pku.edu.cn>
 * @date   Thu Nov 19 12:53:42 2009
 * 
 * @brief  Implementation of neighbor list
 * 
 * 
 */

#define DEVICE_CODE

#include "NeighborList_interface.h"
#include "Auxiliary.h"
#include "NeighborList.h"
#include <stdio.h>
#include "NonBondedInteraction.h"
#include "Reshuffle_interface.h"

/** 
 * these are textures for a fast reference of particle position.
 * 
 */

texture<CoordType, 1, cudaReadModeElementType> global_texRef_neighbor_coord;
texture<TypeType, 1, cudaReadModeElementType> global_texRef_neighbor_type;


void NeighborList::
clearDeviceNeighborList()
{
  if ( mallocedDeviceNeighborList ){
    cudaFree (dnlist.data);
    cudaFree (dnlist.Nneighbor);
    cudaFree (dnlist.forceIndex);
    mallocedDeviceNeighborList = false;
    checkCUDAError ("NeighborList::clearDeviceNeighborList");
  }
}

void NeighborList::
clearNonBondedForce ()
{
  if (mallocedNonBondedForceTable == true){
    cudaFree (nbForceTable);
    mallocedNonBondedForceTable = false;
  }
}

void NeighborList::
clear()
{
  clearDeviceNeighborList();
  clearNonBondedForce();
  unbindGlobalTexture ();
}  

void NeighborList::
unbindGlobalTexture ()
{
  if ( initedGlobalTexture ){
    cudaUnbindTexture(global_texRef_neighbor_coord);
    cudaUnbindTexture(global_texRef_neighbor_type);
    initedGlobalTexture = false;
    checkCUDAError ("NeighborList::unbindGlobalTexture");
  }
}

void NeighborList::
bindGlobalTexture (const MDSystem & sys)
{
  size_t sizetype   = sizeof(TypeType)  *sys.ddata.numMem;
  size_t sizecoord  = sizeof(CoordType) *sys.ddata.numMem;
  cudaBindTexture(0, global_texRef_neighbor_coord, sys.ddata.coord, sizecoord);
  cudaBindTexture(0, global_texRef_neighbor_type, sys.ddata.type, sizetype);
  checkCUDAError ("NeighborList::init texture");
  initedGlobalTexture = true;
}


NeighborList::~NeighborList()
{
  clear();
}

static IndexType hroundUp4 (IndexType x)
{
  if (x & 3 == 0){
    return x;
  }
  else {
    return ((x >> 2) + 1) << 2;
  }
}


void NeighborList::
buildDeviceNeighborListCellList (const MDSystem & sys,
				 const CellList & clist)
{
  dim3 cellBlockDim = clist.getCellBlockDim();
  bool sharednbForceTable (true);
  size_t buildDeviceNeighborList_DeviceCellList_sbuffSize =
      sizeof(IndexType) *	hroundUp4(cellBlockDim.x) +
      sizeof(CoordType) *	hroundUp4(cellBlockDim.x) +
      sizeof(TypeType)  *	hroundUp4(cellBlockDim.x) +
      sizeof(IndexType) *	hroundUp4(nbForceTableLength);
  if (buildDeviceNeighborList_DeviceCellList_sbuffSize >=
      SystemSharedBuffSize - GlobalFunctionParamSizeLimit){
    sharednbForceTable = false;
    buildDeviceNeighborList_DeviceCellList_sbuffSize =
	sizeof(IndexType) *	hroundUp4(cellBlockDim.x) +
	sizeof(CoordType) *	hroundUp4(cellBlockDim.x) +
	sizeof(TypeType)  *	hroundUp4(cellBlockDim.x);
  }
  buildDeviceNeighborList_DeviceCellList 
      <<<clist.getCellGrimDim(), cellBlockDim, 
      buildDeviceNeighborList_DeviceCellList_sbuffSize>>> (
	  // <<<cellGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, 
	  sys.ddata.coord,
	  sys.ddata.type,
	  sys.box,
	  clist.dclist,
	  dnlist,
	  nbForceTable,
	  NatomType,
	  sharednbForceTable,
	  err.ptr_de);
  err.check("NeighborList::buildDeviceNeighborListCellList");
  checkCUDAError ("NeighborList::buildDeviceNeighborListCellList");
}

void NeighborList::
buildDeviceNeighborListAllPair (const MDSystem & sys)
{
  bool sharednbForceTable (true);
  size_t buildDeviceNeighborList_AllPair_sbuffSize =
      sizeof(IndexType) *	hroundUp4(myBlockDim.x) +
      sizeof(CoordType) *	hroundUp4(myBlockDim.x) +
      sizeof(TypeType)  *	hroundUp4(myBlockDim.x) +
      sizeof(IndexType) *	hroundUp4(nbForceTableLength);
  if (buildDeviceNeighborList_AllPair_sbuffSize >=
      SystemSharedBuffSize - GlobalFunctionParamSizeLimit){
    sharednbForceTable = false;
    buildDeviceNeighborList_AllPair_sbuffSize =
	sizeof(IndexType) *	hroundUp4(myBlockDim.x) +
	sizeof(CoordType) *	hroundUp4(myBlockDim.x) +
	sizeof(TypeType)  *	hroundUp4(myBlockDim.x);
  }
  buildDeviceNeighborList_AllPair 
      <<<atomGridDim, myBlockDim,
      buildDeviceNeighborList_AllPair_sbuffSize>>>(
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.type,
	  sys.ddata.rcut,
	  sys.box,
	  dnlist,
	  nbForceTable,
	  NatomType,
	  sharednbForceTable,
	  err.ptr_de);
  err.check("NeighborList::build, build neighbor list all pair");
  checkCUDAError ("NeighborList::build, build neighbor list all pair");
}


void NeighborList::
initNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter)
{
  if (! sysNbInter.beBuilt()) {
    throw MDExcptUnbuiltNonBondedInteraction ("NeighborList");
  }
  NatomType = sysNbInter.numberOfAtomTypes();
  nbForceTableLength = sysNbInter.interactionTableSize();
  cudaMalloc ((void**)&nbForceTable,
	      nbForceTableLength * sizeof(IndexType));
  cudaMemcpy (nbForceTable,
	      sysNbInter.interactionTable(),
	      nbForceTableLength * sizeof(IndexType),
	      cudaMemcpyHostToDevice);
  checkCUDAError ("AtomNBForceTable::deviceInitTable");
  mallocedNonBondedForceTable = true;
}

void NeighborList::
mallocDeviceNeighborList (const MDSystem & sys,
			  const ScalorType & DeviceNeighborListExpansion)
{
  ScalorType density = sys.ddata.numAtom / (sys.box.size.x * sys.box.size.y * sys.box.size.z);
  ScalorType expectedNumberInList 
      = 4./3. * M_PI * myrlist * myrlist * myrlist * density;
  dnlist.listLength = IndexType(expectedNumberInList * DeviceNeighborListExpansion);
  if (dnlist.listLength < 30){
    dnlist.listLength = 30;
  }
  printf ("#@ length of the neighbor list is %d\n", dnlist.listLength);
  cudaMalloc ((void**)&(dnlist.data), sizeof(IndexType) * dnlist.stride * dnlist.listLength);
  cudaMalloc ((void**)&(dnlist.Nneighbor), sizeof(IndexType) * sys.ddata.numAtom);
  cudaMalloc ((void**)&(dnlist.forceIndex), sizeof(IndexType) *  dnlist.stride * dnlist.listLength);
  // reshuffle backup things
  cudaMalloc ((void**)&(bkdnlistData), sizeof(IndexType) * dnlist.stride * dnlist.listLength);
  cudaMalloc ((void**)&(bkdnlistNneighbor), sizeof(IndexType) * sys.ddata.numAtom);
  cudaMalloc ((void**)&(bkdnlistForceIndex), sizeof(IndexType) *  dnlist.stride * dnlist.listLength);
  checkCUDAError ("NeighborList::mallocDeviceNeighborList");
  mallocedDeviceNeighborList = true;
}



void NeighborList::
reinit (const SystemNonBondedInteraction & sysNbInter,
	const MDSystem & sys,
	const ScalorType & rlist,
	const ScalorType & rlistExten,
	const IndexType & NTread,
	const ScalorType & DeviceNeighborListExpansion)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;

  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);
  
  myrlist = rlist;
  dnlist.rlist = myrlist;
  dnlist.rlistExten = rlistExten;
  dnlist.stride = sys.ddata.numAtom;
  
  // init neighbor list
  clearDeviceNeighborList ();
  mallocDeviceNeighborList (sys, DeviceNeighborListExpansion);  

  clearNonBondedForce ();
  initNonBondedInteraction (sysNbInter);

  unbindGlobalTexture ();
  bindGlobalTexture (sys);
  
  //init shared memory size
}

void NeighborList::
rebuild (const MDSystem & sys,
	 const CellList & clist,
	 MDTimer * timer)
{
  if (clist.isempty()){
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    // printf ("rlist is %f\n", dnlist.rlist);
    buildDeviceNeighborListAllPair (sys);
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }
  else {
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    buildDeviceNeighborListCellList (sys, clist);
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }  
}


void NeighborList::
reshuffle (const IndexType * indexTable,
	   const IndexType & numAtom,
	   MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeReshuffleSystem);

  // Reshuffle_reshuffleDeviceCellList
  //     <<<cellGridDim, myBlockDim>>> (
  // 	  dclist.data, indexTable);
  // cudaMemcpy (bkbackupCoord, backupCoord,
  // 	      sizeof (CoordType) * numAtom,
  // 	      cudaMemcpyDeviceToDevice);
  // Reshuffle_reshuffleArray
  //     <<<atomGridDim, myBlockDim>>> 
  //     (bkbackupCoord, numAtom, indexTable, backupCoord);
  Reshuffle_backupDeviceNeighborList
      <<<atomGridDim, myBlockDim,
      2 * myBlockDim.x * sizeof(IndexType)>>> (
	  numAtom,
	  dnlist.data,
	  dnlist.forceIndex,
	  dnlist.stride,
	  dnlist.Nneighbor,
	  bkdnlistData,
	  bkdnlistForceIndex,
	  bkdnlistNneighbor);
  checkCUDAError ("NeighborList::reshuffle backup");
  Reshuffle_reshuffleDeviceNeighborList
      <<<atomGridDim, myBlockDim,
      2 * myBlockDim.x * sizeof(IndexType)>>> (
	  numAtom,
	  bkdnlistData,
	  bkdnlistForceIndex,
	  dnlist.stride,
	  bkdnlistNneighbor,
	  indexTable,
	  dnlist.data,
	  dnlist.forceIndex,
	  dnlist.Nneighbor);
  checkCUDAError ("NeighborList::reshuffle reshuffle");
  if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
}

NeighborList::
NeighborList (const SystemNonBondedInteraction & sysNbInter,
	      const MDSystem & sys,
	      const ScalorType & rlist,
	      const ScalorType & rlistExten,
	      const IndexType & NTread,
	      const ScalorType & DeviceNeighborListExpansion)
    : mallocedDeviceNeighborList (false),
      mallocedNonBondedForceTable (false),
      initedGlobalTexture (false)
{
  reinit (sysNbInter, sys, rlist, rlistExten, NTread, DeviceNeighborListExpansion);
}





////////////////////////////////////////////////////////////
// for the reason of using texture, we place this function here. it
// should be placed in NeighborList.cu
////////////////////////////////////////////////////////////

using namespace RectangularBoxGeometry;

__device__ IndexType
shiftedD3toD1 (DeviceCellList clist,
	       RectangularBox box,
	       int ix,
	       int iy,
	       int iz,
	       ScalorType * shiftx ,
	       ScalorType * shifty,
	       ScalorType * shiftz)
{
  int tmp;
  ix += (tmp = -int(floorf(ix * clist.NCelli.x))) * clist.NCell.x;
  *shiftx = tmp * box.size.x;
  iy += (tmp = -int(floorf(iy * clist.NCelli.y))) * clist.NCell.y;
  *shifty = tmp * box.size.y;
  iz += (tmp = -int(floorf(iz * clist.NCelli.z))) * clist.NCell.z;
  *shiftz = tmp * box.size.z;
  return D3toD1 (clist.NCell, ix, iy, iz);
}

__global__ void
buildDeviceNeighborList_DeviceCellList (const IndexType		numAtom,
					const CoordType *	coord,
					const TypeType *	type,
					const RectangularBox	box,
					const DeviceCellList	clist,
					DeviceNeighborList	nlist,
					const IndexType *	nbForceTable,
					const IndexType		NatomType,
					const bool		sharednbForceTable,
					mdError_t *		ptr_de )
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
  // set number of neighbor to 0
  IndexType Nneighbor = 0;
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_neighbor_coord, ii);
    reftype = tex1Dfetch(global_texRef_neighbor_type, ii);
#endif
  }
  ScalorType rlist = nlist.rlist;

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];
  IndexType * nbForceTableBuff = NULL;

  IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  if (sharednbForceTable){
    nbForceTableBuff = (IndexType *) &targettype[roundUp4(blockDim.x)];
    cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  }
  __syncthreads();

  // __shared__ volatile  IndexType  targetIndexes [MaxThreadsPerBlock];
  // __shared__ volatile CoordType target [MaxThreadsPerBlock];
  // __shared__ volatile  TypeType   targettype    [MaxThreadsPerBlock];
  // __shared__ volatile  IndexType nbForceTableBuff [MaxNBForceTableBuffSize];

  // IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  // if (sharednbForceTable){
  //   cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  // }
  // __syncthreads();

  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (clist.NCell.x == 1) oneCellX = true;
  if (clist.NCell.y == 1) oneCellY = true;
  if (clist.NCell.z == 1) oneCellZ = true;
  ScalorType rlist2 = rlist * rlist;
  
  // loop over 27 neighbor cells
  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    // if (threadIdx.x == 0){
    //   printf ("%d %d\n", bid, clist.numNeighborCell[bid]);
    // }
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex    (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_neighbor_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_neighbor_type, targetIndexes[tid]);
    }
    __syncthreads();
	
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
	    targetIndexes[jj] != ii){
	  IndexType fidx;
	  if (sharednbForceTable){
	    fidx = AtomNBForceTable::calForceIndex (
		nbForceTableBuff, NatomType, reftype, targettype[jj]);
	  }
	  else {
	    fidx = AtomNBForceTable::calForceIndex (
		nbForceTable, NatomType, reftype, targettype[jj]);
	  }
	  // if (fidx != mdForceNULL) {
	  IndexType listIdx = Nneighbor * nlist.stride + ii;
	  nlist.data[listIdx] = targetIndexes[jj];
	  nlist.forceIndex[listIdx] = fidx;
	  Nneighbor ++;
	  // }
	}
      }
    }
  }

  if (ii != MaxIndexValue) {
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
      return;
    }
    nlist.Nneighbor[ii] = Nneighbor;
    // printf ("%d %d\n", ii, Nneighbor);
  }
}



__global__ void
buildDeviceCellList_step1 (IndexType		numAtom,
			   CoordType *		coord,
			   IntScalorType *	coordNoix,
			   IntScalorType *	coordNoiy,
			   IntScalorType *	coordNoiz,
			   RectangularBox	box,
			   DeviceCellList	clist,
			   IndexType *		sendBuff,
			   IndexType *		targetBuff,
			   mdError_t *		ptr_de,
			   IndexType *		erridx,
			   ScalorType *		errsrc)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  extern __shared__ volatile IndexType sbuff[];
  volatile IndexType * originalData = (volatile IndexType *) sbuff;
  volatile IndexType * targetCellid = (volatile IndexType *) &originalData[blockDim.x];
  
  // __shared__ volatile  IndexType originalData[MaxThreadsPerBlock];
  // __shared__ volatile  IndexType targetCellid[MaxThreadsPerBlock];
  
  // copy data from cell list
  originalData[tid] = clist.data[bid*clist.stride + tid];
  IndexType originalNumber = clist.numbers[bid];

  // calculate the target cell
  if (originalData[tid] != MaxIndexValue){
    IndexType targetCelli, targetCellj, targetCellk;
    IndexType thisid = originalData[tid];
#ifdef COMPILE_NO_TEX
    ref = coord[thisid];
#else
    CoordType ref (tex1Dfetch(global_texRef_neighbor_coord, thisid));
#endif
    targetCelli = IndexType(ref.x * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(ref.y * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(ref.z * box.sizei.z * ScalorType (clist.NCell.z));
    if (targetCelli == clist.NCell.x){
      targetCelli -= clist.NCell.x;
      coord[thisid].x -= box.size.x;
      coordNoix[thisid] ++;
    }
    if (targetCellj == clist.NCell.y){
      targetCellj -= clist.NCell.y;
      coord[thisid].y -= box.size.y;
      coordNoiy[thisid] ++;
    }
    if (targetCellk == clist.NCell.z){
      targetCellk -= clist.NCell.z;
      coord[thisid].z -= box.size.z;
      coordNoiz[thisid] ++;
    }
    targetCellid[tid] = D3toD1 (clist.NCell, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = ref.x;
	// return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = ref.y;
	// return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = ref.z;
	// return;
      }      
    }
  }
  else {
    targetCellid[tid] = MaxIndexValue;
  }

  // mark particles to be send 
  IndexType mark = MaxIndexValue - (MaxIndexValue >> 1);
  if (tid < originalNumber && targetCellid[tid] != bid){
    originalData[tid] += mark;
  }
  
  // head sort
  IndexType total1 = headSort (originalData, targetCellid);
  IndexType total0 = blockDim.x - total1;
  
  // unmark and copy to send buff
  if (tid < originalNumber && targetCellid[tid] != bid){
    sendBuff  [bid*clist.stride + tid - total0] = originalData[tid] - mark;
    targetBuff[bid*clist.stride + tid - total0] = targetCellid[tid];
    originalData[tid] = MaxIndexValue;
  }
  __syncthreads();
  
  // modify cell list
  clist.data[bid*clist.stride + tid] = originalData[tid];
  if (tid == 0) clist.numbers[bid] = total0;
}



__global__ void
buildDeviceCellList_initBuff (IndexType * sendBuff,
			      IndexType * targetBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
  targetBuff[ii] = 0;
}

__global__ void
buildDeviceCellList_clearBuff (IndexType * sendBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
}

  
  
__global__ void
buildDeviceCellList_step2 (RectangularBox	box,
			   DeviceCellList	clist,
			   IndexType *		sendBuff,
			   IndexType *		targetBuff,
			   IndexType		bitDeepth,
			   mdError_t *		ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType thisid;
  IndexType ii = 0;
  IndexType buffPosi;
  if (tid == 0){
    while ((thisid = sendBuff[buffPosi = (bid*clist.stride + ii)]) != MaxIndexValue){
      IndexType cellid = targetBuff[buffPosi];
      IndexType tailIdx = atomicInc(&clist.numbers[cellid], blockDim.x);
      if (tailIdx >= blockDim.x&& ptr_de != NULL) {
	*ptr_de = mdErrorShortCellList;
	return;
      } 
      clist.data[cellid * clist.stride + tailIdx] = thisid;
      sendBuff[buffPosi] = MaxIndexValue;
      ii ++;
    }
  }
}






