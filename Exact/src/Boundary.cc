#include "cctk.h"
#include "cctk_Arguments.h"
#include "PreSync.h"

void Bndry_Exact_extrap(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                 CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *table_handles);

void Exact_RegisterBC(CCTK_ARGUMENTS)
{

  int err = 0;
  err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_Exact_extrap, "Exact_linear_extrap");
  if (err) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error %d when registering routine to handle \"Scalar\" "
               "boundary condition",
               err);
  }
}

void Exact_SelectBCs(CCTK_ARGUMENTS)
{
  int ierr = 0;

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "ADMBase::metric", "Exact_linear_extrap");
  if (ierr < 0)
     CCTK_ERROR("Failed to register Exact_linear_extrap BC for ADMBase::metric!");

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "ADMBase::curv", "Exact_linear_extrap");
  if (ierr < 0)
     CCTK_ERROR("Failed to register Exact_linear_extrap BC for ADMBase::curv!");

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "Exact::Exact_slice", "Exact_linear_extrap");
  if (ierr < 0)
     CCTK_ERROR("Failed to register Exact_linear_extrap BC for Exact::Exact_slice!");

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "Exact::Exact_slicetemp1", "Exact_linear_extrap");
  if (ierr < 0)
     CCTK_ERROR("Failed to register Exact_linear_extrap BC for Exact::Exact_slicetemp1!");

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "Exact::Exact_slicetemp2", "Exact_linear_extrap");
  if (ierr < 0)
     CCTK_ERROR("Failed to register Exact_linear_extrap BC for Exact::Exact_slicetemp2!");
}

void Bndry_Exact_extrap(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                 CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *table_handles)
{
  int i,j,k;
  int nx = cctkGH->cctk_lsh[0]-1;
  int ny = cctkGH->cctk_lsh[1]-1;
  int nz = cctkGH->cctk_lsh[2]-1;

  CCTK_REAL *var = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, CCTK_FullName(var_indices[0]));

  //Linear extrapolation from the interiors to the boundaries.
  //Does not support octant or quadrant.

  //6 faces.
  if (cctkGH->cctk_bbox[1] == 1 && nx >= 4 ) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
        var[0,j,k] = 2 * var[1,j,k] - var[2,j,k];
      }
    }
  }
  if (cctkGH->cctk_bbox[2] == 1 && nx >= 4 ) {
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
        var[nx,j,k] = 2 * var[nx-1,j,k] - var[nx-2,j,k];
      }
    }
  }
  if (cctkGH->cctk_bbox[3] == 1 && ny >= 4 ) {
    for(i = 0; i < nx; i++){
      for(k = 0; k < nz; k++){
        var[i,0,k] = 2 * var[i,1,k] - var[i,2,k];
      }
    }
  }
  if (cctkGH->cctk_bbox[4] == 1 && ny >= 4 ) {
    for(i = 0; i < nx; i++){
      for(k = 0; k < nz; k++){
        var[i,ny,k] = 2 * var[i,ny-1,k] - var[i,ny-2,k];
      }
    }
  }
  if (cctkGH->cctk_bbox[5] == 1 && nz >= 4 ) {
    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
        var[i,j,0] = 2 * var[i,j,1] - var[i,j,2];
      }
    }
  }
  if (cctkGH->cctk_bbox[6] == 1 && nz >= 4 ) {
    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
        var[i,j,nz] = 2 * var[i,j,nz-1] - var[i,j,nz-2];
      }
    }
  }

  //12 edges.
  //4 around face x=min.
  if (cctkGH->cctk_bbox[1] == 1 && nx >= 4 && cctkGH->cctk_bbox[3] == 1 && ny >= 4) {
    for(k = 0; k < nz; k++){
      var[0,0,k] = 2 * var[1,1,k] - var[2,2,k];
    }
  }
  if (cctkGH->cctk_bbox[1] == 1 && nx >= 4 && cctkGH->cctk_bbox[4] == 1 && ny >= 4) {
    for(k = 0; k < nz; k++){
      var[0,ny,k] = 2 * var[1,ny-1,k] - var[2,ny-2,k];
    }
  }
  if (cctkGH->cctk_bbox[1] == 1 && nx >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(j = 0; j < ny; j++){
      var[0,j,0] = 2 * var[1,j,1] - var[2,j,2];
    }
  }
  if (cctkGH->cctk_bbox[1] == 1 && nx >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(j = 0; j < ny; j++){
      var[0,j,nz] = 2 * var[1,j,nz-1] - var[2,j,nz-2];
    }
  }

  //4 around face x=max.
  if (cctkGH->cctk_bbox[2] == 1 && nx >= 4 && cctkGH->cctk_bbox[3] == 1 && ny >= 4) {
    for(k = 0; k < nz; k++){
      var[nx,0,k] = 2 * var[nx-1,1,k] - var[nx-2,2,k];
    }
  }
  if (cctkGH->cctk_bbox[2] == 1 && nx >= 4 && cctkGH->cctk_bbox[4] == 1 && ny >= 4) {
    for(k = 0; k < nz; k++){
      var[nx,ny,k] = 2 * var[nx-1,ny-1,k] - var[nx-2,ny-2,k];
    }
  }
  if (cctkGH->cctk_bbox[2] == 1 && nx >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(j = 0; j < ny; j++){
      var[nx,j,0] = 2 * var[nx-1,j,1] - var[nx-2,j,2];
    }
  }
  if (cctkGH->cctk_bbox[2] == 1 && nx >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(j = 0; j < ny; j++){
      var[nx,j,nz] = 2 * var[nx-1,j,nz-1] - var[nx-2,j,nz-2];
    }
  }

  //Remaining 2 in y=min.
  if (cctkGH->cctk_bbox[3] == 1 && ny >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(i = 0; i < nx; i++){
      var[i,0,0] = 2 * var[i,1,1] - var[i,2,2];
    }
  }
  if (cctkGH->cctk_bbox[3] == 1 && ny >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(i = 0; i < nx; i++){
      var[i,0,nz] = 2 * var[i,1,nz-1] - var[i,2,nz-2];
    }
  }

  //Remaining 2 in y=ymax.      
  if (cctkGH->cctk_bbox[4] == 1 && ny >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(i = 0; i < nx; i++){
      var[i,ny,0] = 2 * var[i,ny-1,1] - var[i,ny-2,2];
    }
  }
  if (cctkGH->cctk_bbox[4] == 1 && ny >= 4 && cctkGH->cctk_bbox[5] == 1 && nz >= 4) {
    for(i = 0; i < nx; i++){
      var[i,ny,nz] = 2 * var[i,ny-1,nz-1] - var[i,ny-2,nz-2];
    }
  }

  //8 corners.
  if(nx >= 4 && ny >= 4 && nz >= 4) {
    if (cctkGH->cctk_bbox[1] == 1 && cctkGH->cctk_bbox[3] == 1 && cctkGH->cctk_bbox[5] == 1) {
      var[0,0,0] = 2 * var[1,1,1] - var[2,2,2];
    }
    if (cctkGH->cctk_bbox[1] == 1 && cctkGH->cctk_bbox[3] == 1 && cctkGH->cctk_bbox[6] == 1) {
      var[0,0,nz] = 2 * var[1,1,nz-1] - var[2,2,nz-2];
    }
    if (cctkGH->cctk_bbox[1] == 1 && cctkGH->cctk_bbox[4] == 1 && cctkGH->cctk_bbox[5] == 1) {
      var[0,ny,0] = 2 * var[1,ny-1,1] - var[2,ny-2,2];
    }
    if (cctkGH->cctk_bbox[1] == 1 && cctkGH->cctk_bbox[4] == 1 && cctkGH->cctk_bbox[6] == 1) {
      var[0,ny,nz] = 2 * var[1,ny-1,nz-1] - var[2,ny-2,nz-2];
    }
    if (cctkGH->cctk_bbox[2] == 1 && cctkGH->cctk_bbox[3] == 1 && cctkGH->cctk_bbox[5] == 1) {
      var[nx,0,0] = 2 * var[nx-1,1,1] - var[nx-2,2,2];
    }
    if (cctkGH->cctk_bbox[2] == 1 && cctkGH->cctk_bbox[3] == 1 && cctkGH->cctk_bbox[6] == 1) {
      var[nx,0,nz] = 2 * var[nx-1,1,nz-1] - var[nx-2,2,nz-2];
    }
    if (cctkGH->cctk_bbox[2] == 1 && cctkGH->cctk_bbox[4] == 1 && cctkGH->cctk_bbox[5] == 1) {
      var[nx,ny,0] = 2 * var[nx-1,ny-1,1] - var[nx-2,ny-2,2];
    }
    if (cctkGH->cctk_bbox[2] == 1 && cctkGH->cctk_bbox[4] == 1 && cctkGH->cctk_bbox[6] == 1) {
      var[nx,ny,nz] = 2 * var[nx-1,ny-1,nz-1] - var[nx-2,ny-2,nz-2];
    }
  }
}
