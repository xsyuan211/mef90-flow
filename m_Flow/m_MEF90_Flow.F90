#include "../MEF90/mef90.inc"
Module m_MEF90_Flow
#include "finclude/petscdef.h"
#include "finclude/petsctsdef.h"
   Use m_MEF90
   Use m_MEF90_FlowCtx
   Use m_MEF90_DefMechCtx_Type   !!! for MEF90DefMechCtx_Type reference
   Use m_MEF90_FlowAssembly2D, &   !!! NEW: bridge to 2D dimmed assembly template
      MEF90FlowOperator2D               => MEF90FlowOperator, &
      MEF90FlowBilinearForm2D           => MEF90FlowBilinearForm, &
      MEF90FlowIFunction2D              => MEF90FlowIFunction, &
      MEF90FlowIJacobian2D              => MEF90FlowIJacobian, &
      MEF90FlowEnergy2D                 => MEF90FlowEnergy, &
      MEF90FlowViewEXO2D                => MEF90FlowViewEXO, &
      MEF90FlowUpdateboundaryPressure2D => MEF90FlowUpdateboundaryPressure, &
      MEF90FlowSetTransients2D          => MEF90FlowSetTransients, &
      MEF90FlowComputePressureGradient2D => MEF90FlowComputePressureGradient, &
      MEF90FlowApplyBiotForce2D         => MEF90FlowApplyBiotForce
   Use m_MEF90_FlowAssembly3D, &   !!! NEW: bridge to 3D dimmed assembly template
      MEF90FlowOperator3D               => MEF90FlowOperator, &
      MEF90FlowBilinearForm3D           => MEF90FlowBilinearForm, &
      MEF90FlowIFunction3D              => MEF90FlowIFunction, &
      MEF90FlowIJacobian3D              => MEF90FlowIJacobian, &
      MEF90FlowEnergy3D                 => MEF90FlowEnergy, &
      MEF90FlowViewEXO3D                => MEF90FlowViewEXO, &
      MEF90FlowUpdateboundaryPressure3D => MEF90FlowUpdateboundaryPressure, &
      MEF90FlowSetTransients3D          => MEF90FlowSetTransients, &
      MEF90FlowComputePressureGradient3D => MEF90FlowComputePressureGradient, &
      MEF90FlowApplyBiotForce3D         => MEF90FlowApplyBiotForce

   Implicit none

   Private
   !!! Public — Context lifecycle (defined in FlowCtx.F90 or here)
   Public :: MEF90FlowCtxSetFromOptions,     &
             MEF90FlowCtxCreateVectors,      &
             MEF90FlowCtxSetSections,        &
             MEF90FlowCtxDestroyVectors,     &
             MEF90FlowCtxDestroy,            &
             MEF90FlowUpdatePerm,            &   !!! k(v) = k_min + (k_max - k_min)*damage
             MEF90FlowUpdateStorativity,     &   !!! S per cell from options
             MEF90FlowCreateSNES,            &   !!! SNES nonlinear pressure solver
             MEF90FlowCreateTS,              &   !!! TS transient pressure solver
             MEF90FlowSetTransients,         &   !!! time-dependent source/BC update
             MEF90FlowUpdateboundaryPressure,&
             MEF90FlowOperator,              &   !!! Assembled via 2D/3D dispatch below
             MEF90FlowBilinearForm,          &
             MEF90FlowIFunction,             &
             MEF90FlowIJacobian,             &
             MEF90FlowEnergy,                &
             MEF90FlowViewEXO,               &
             MEF90FlowComputePressureGradient,& !!! Cell-centered ∇p volume average
             MEF90FlowApplyBiotForce           !!! Biot f_p = -α·∇p → DefMechCtx%force

!!! --------------------------------------------------------------------------------
!!!  Wrapper functions — dispatch 2D vs 3D based on mesh dimension (HeatXfer pattern)
!!! --------------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowOperator"
   Subroutine MEF90FlowOperator(snesPressure,x,residual,FlowCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesPressure
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowOperator2D(snesPressure,x,residual,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowOperator3D(snesPressure,x,residual,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowOperator

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowBilinearForm"
   Subroutine MEF90FlowBilinearForm(snesPressure,x,Jac,PreJac,residualMode,FlowCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesPressure
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      PetscReal,Intent(IN)                               :: residualMode
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowBilinearForm2D(snesPressure,x,Jac,PreJac,residualMode,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowBilinearForm3D(snesPressure,x,Jac,PreJac,residualMode,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIFunction"
   Subroutine MEF90FlowIFunction(ts,t,x,xDot,F,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: ts
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      Type(Vec),Intent(INOUT)                            :: F
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowIFunction2D(ts,t,x,xDot,F,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowIFunction3D(ts,t,x,xDot,F,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowIFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowIJacobian"
   Subroutine MEF90FlowIJacobian(ts,t,x,xDot,a,Jac,PreJac,FlowCtx,ierr)
      Type(TS),Intent(IN)                                :: ts
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x, xDot
      PetscReal,Intent(IN)                               :: a
      Type(Mat),Intent(INOUT)                            :: Jac, PreJac
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowIJacobian2D(ts,t,x,xDot,a,Jac,PreJac,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowIJacobian3D(ts,t,x,xDot,a,Jac,PreJac,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowIJacobian

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowEnergy"
   Subroutine MEF90FlowEnergy(pressure,timeVal,FlowCtx,Energy,InjectionWork,ierr)
      Type(Vec),Intent(IN)                               :: pressure
      PetscReal,Intent(IN)                               :: timeVal
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscReal,Intent(OUT)                              :: Energy, InjectionWork
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowEnergy2D(pressure,timeVal,FlowCtx,Energy,InjectionWork,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowEnergy3D(pressure,timeVal,FlowCtx,Energy,InjectionWork,ierr)
      End If
   End Subroutine MEF90FlowEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowViewEXO"
   Subroutine MEF90FlowViewEXO(FlowCtx,step,ierr)
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowViewEXO2D(FlowCtx,step,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowViewEXO3D(FlowCtx,step,ierr)
      End If
   End Subroutine MEF90FlowViewEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdateboundaryPressure"
   Subroutine MEF90FlowUpdateboundaryPressure(boundaryPressure,FlowCtx,ierr)
      Type(Vec),Intent(INOUT)                            :: boundaryPressure
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowUpdateboundaryPressure2D(boundaryPressure,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowUpdateboundaryPressure3D(boundaryPressure,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowUpdateboundaryPressure

#undef __FUNCT__
#define __FUNCT__ "MEF90FlowSetTransients"
   Subroutine MEF90FlowSetTransients(FlowCtx,step,timeVal,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      PetscInt,Intent(IN)                                :: step
      PetscReal,Intent(IN)                               :: timeVal
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowSetTransients2D(FlowCtx,step,timeVal,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowSetTransients3D(FlowCtx,step,timeVal,ierr)
      End If
   End Subroutine MEF90FlowSetTransients

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Update permeability field from damage field and global options
!!!  k(v) = k_min + (k_max - k_min) * max(v, residual_damage)
!!!  permField(cell) = k(v) / mu   (cell-centered, averaged from damage)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdatePerm"

   Subroutine MEF90FlowUpdatePerm(FlowCtx,DefMechCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: DefMechCtx
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellOpt
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer     :: elem
      Type(MEF90Element_Type)                             :: elemType
      Type(IS)                                            :: CellSetGlobalIS, setIS, setIdx
      PetscInt,dimension(:),Pointer                       :: setID
      PetscInt                                            :: set, cell
      PetscReal                                           :: damageGauss, vGauss, k_v, k_eff
      PetscReal,Dimension(:),Pointer                      :: dmgPtr
      PetscReal,Dimension(:),Pointer                      :: permPtr
      PetscInt                                            :: iDoF, iGauss, QuadratureOrder

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellOpt,ierr);CHKERRQ(ierr)

         !!! Base k for this cell set — use override if given, else k_min
         k_eff = gOpt%permeabilityMin
         If (cellOpt%permeabilityOverride > 0.0_Kr) Then
            k_eff = cellOpt%permeabilityOverride
         End If

         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)

         elemType = MEF90knownElements(cellOpt%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(FlowCtx%DM,setIS,elem,QuadratureOrder, &
              cellOpt%ElemTypeShortID,ierr);CHKERRQ(ierr)

         Allocate(dmgPtr(elemType%numDof))
         Do cell = 1, size(setIdx)
            !!! Restrict damage to this cell's DOFs
            Call SectionRealRestrictClosure(DefMechCtx%DMScalSec,FlowCtx%DM,setIdx(cell),elemType%numDof,dmgPtr,ierr);CHKERRQ(ierr)

            !!! Average damage at cell (Gauss-point averaging of v-field shape)
            vGauss = 0.0_Kr
            Do iGauss = 1, size(elem(cell)%Gauss_C)
               damageGauss = 0.0_Kr
               Do iDoF = 1, elemType%numDof
                  damageGauss = damageGauss + dmgPtr(iDoF) * elem(cell)%BF(iDoF,iGauss)
               End Do
               vGauss = vGauss + elem(cell)%Gauss_C(iGauss) * damageGauss
            End Do
            vGauss = vGauss / elem(cell)%volume

            !!! Apply residual damage floor
            If (vGauss < gOpt%residualDamage) vGauss = gOpt%residualDamage

            !!! k(v) = k_min + (k_max - k_min)*v
            k_v = k_eff + (gOpt%permeabilityMax - k_eff) * vGauss

            !!! Store as cell scalar: permField(cell) = k_v / mu
            Allocate(permPtr(1))
            Call SectionRealRestrict(FlowCtx%cellDMScalSec,setIdx(cell),permPtr,ierr);CHKERRQ(ierr)
            permPtr(1) = k_v / max(gOpt%dynamicViscosity,1e-12_Kr)
            Call SectionRealRestore(FlowCtx%cellDMScalSec,setIdx(cell),permPtr,ierr);CHKERRQ(ierr)
            DeAllocate(permPtr)
         End Do
         DeAllocate(dmgPtr)
         Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do

      Call SectionRealComplete(FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(FlowCtx%cellDMScalSec,FlowCtx%cellDMScalScatter, &
           SCATTER_FORWARD,FlowCtx%permField,ierr);CHKERRQ(ierr)

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowUpdatePerm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Storativity per cell — store into FlowCtx%storativity Vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowUpdateStorativity"

   Subroutine MEF90FlowUpdateStorativity(FlowCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt
      Type(MEF90FlowCellSetOptions_Type),pointer         :: cellOpt
      Type(IS)                                            :: CellSetGlobalIS, setIS, setIdx
      PetscInt,dimension(:),Pointer                       :: setID
      PetscInt                                            :: set, cell
      PetscReal                                           :: S_eff
      PetscReal,Dimension(:),Pointer                      :: Sptr

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(FlowCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1, size(setID)
         Call PetscBagGetDataMEF90FlowCtxCellSetOptions(FlowCtx%CellSetOptionsBag(set),cellOpt,ierr);CHKERRQ(ierr)
         S_eff = gOpt%storativity
         If (cellOpt%storativityOverride >= 0.0_Kr) S_eff = cellOpt%storativityOverride

         Call DMMeshGetStratumIS(FlowCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Do cell = 1, size(setIdx)
            Allocate(Sptr(1))
            Call SectionRealRestrict(FlowCtx%cellDMScalSec,setIdx(cell),Sptr,ierr);CHKERRQ(ierr)
            Sptr(1) = S_eff
            Call SectionRealRestore(FlowCtx%cellDMScalSec,setIdx(cell),Sptr,ierr);CHKERRQ(ierr)
            DeAllocate(Sptr)
         End Do
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do

      Call SectionRealComplete(FlowCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(FlowCtx%cellDMScalSec,FlowCtx%cellDMScalScatter, &
           SCATTER_FORWARD,FlowCtx%storativity,ierr);CHKERRQ(ierr)

      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowUpdateStorativity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CREATE SNES NONLINEAR SOLVER for pressure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCreateSNES"

   Subroutine MEF90FlowCreateSNES(FlowCtx,SNESPressure,residual,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(SNES),Intent(OUT)                              :: SNESPressure
      Type(Vec),Intent(INOUT)                             :: residual
      PetscErrorCode,Intent(OUT)                          :: ierr

      Type(MEF90FlowGlobalOptions_Type),pointer          :: gOpt

      Call PetscBagGetDataMEF90FlowCtxGlobalOptions(FlowCtx%GlobalOptionsBag,gOpt,ierr);CHKERRQ(ierr)

      Call SNESCreate(FlowCtx%MEF90Ctx%comm,SNESPressure,ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(SNESPressure,ierr);CHKERRQ(ierr)
      Call SNESSetDM(SNESPressure,FlowCtx%DMScal,ierr);CHKERRQ(ierr)
      Call SNESSetJacobian(SNESPressure,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
           MEF90FlowBilinearForm,FlowCtx,ierr);CHKERRQ(ierr)
      Call SNESSetFunction(SNESPressure,residual,MEF90FlowOperator,FlowCtx,ierr);CHKERRQ(ierr)
      Call SNESSetConvergence(SNESPressure,gOpt%pressureRTol,gOpt%pressureATol, &
           PETSC_DEFAULT_REAL,gOpt%pressureMaxIt,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(SNESPressure,"SNESPressure",ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCreateSNES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  CREATE TS TRANSIENT SOLVER for pressure (theta-method/BDF)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowCreateTS"

   Subroutine MEF90FlowCreateTS(FlowCtx,tsPressure,residual,timeStart,dt,ierr)
      Type(MEF90FlowCtx_Type),Intent(INOUT)              :: FlowCtx
      Type(TS),Intent(OUT)                                :: tsPressure
      Type(Vec),Intent(INOUT)                             :: residual
      PetscReal,Intent(IN)                                :: timeStart, dt
      PetscErrorCode,Intent(OUT)                          :: ierr

      Call TSCreate(FlowCtx%MEF90Ctx%comm,tsPressure,ierr);CHKERRQ(ierr)
      Call TSSetFromOptions(tsPressure,ierr);CHKERRQ(ierr)
      Call TSSetType(tsPressure,TSBEULER,ierr);CHKERRQ(ierr)
      Call TSSetDM(tsPressure,FlowCtx%DMScal,ierr);CHKERRQ(ierr)
      Call TSSetIFunction(tsPressure,residual,MEF90FlowIFunction,FlowCtx,ierr);CHKERRQ(ierr)
      Call TSSetIJacobian(tsPressure,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
           MEF90FlowIJacobian,FlowCtx,ierr);CHKERRQ(ierr)
      Call TSSetTime(tsPressure,timeStart,ierr);CHKERRQ(ierr)
      Call TSSetTimeStep(tsPressure,dt,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(tsPressure,"TSPressure",ierr);CHKERRQ(ierr)
   End Subroutine MEF90FlowCreateTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  COMPUTE CELL-CENTERED ∇p — volume-weighted average per cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowComputePressureGradient"

   Subroutine MEF90FlowComputePressureGradient(pressure,cellPressureGradient,FlowCtx,ierr)
      Type(Vec),Intent(IN)                               :: pressure
      Type(Vec),Intent(INOUT)                            :: cellPressureGradient
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscInt                                           :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowComputePressureGradient2D(pressure,cellPressureGradient,FlowCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowComputePressureGradient3D(pressure,cellPressureGradient,FlowCtx,ierr)
      End If
   End Subroutine MEF90FlowComputePressureGradient


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  APPLY BIOT BODY FORCE  f_p = -α · ∇p  →  DefMechCtx%force
!!!  Internal dim dispatch — ComputePressureGradient depends on dim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90FlowApplyBiotForce"

   Subroutine MEF90FlowApplyBiotForce(FlowCtx,DefMechCtx,ierr)
      Type(MEF90FlowCtx_Type),Intent(IN)                 :: FlowCtx
      Type(MEF90DefMechCtx_Type),Intent(INOUT)           :: DefMechCtx
      PetscErrorCode,Intent(OUT)                          :: ierr
      PetscInt                                            :: dim
      Call DMMeshGetDimension(FlowCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90FlowApplyBiotForce2D(FlowCtx,DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90FlowApplyBiotForce3D(FlowCtx,DefMechCtx,ierr)
      End If
   End Subroutine MEF90FlowApplyBiotForce

End Module m_MEF90_Flow
