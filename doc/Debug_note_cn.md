---
title: ATEK的移植
---

[TOC]

# AETK简介

AETK是HIT的Chun Shen开发的求解器。

# 移植



从AETK的说明来看，除了安装OpenFOAM 2.0.1，它还需要重新编译`finiteVolume`库，其中涉及到替换里面的`MRF`库。

上网查看了一下[OpenFOAM 2.0.x的代码](https://github.com/OpenFOAM/OpenFOAM-2.0.x)，其中的

## AETKv1.0的MRF库

### 文件

对比三组文件：

- [AETKv1/src/timeStepping/MRF/](AETKv1/src/timeStepping/MRF/)
- [OpenFOAM-2.0.x/src/finiteVolume/cfdTools/general/MRF/](OpenFOAM-2.0.x/src/finiteVolume/cfdTools/general/MRF/)
- [OpenFOAM-4.x/src/finiteVolume/cfdTools/general/MRF/](OpenFOAM-4.x/src/finiteVolume/cfdTools/general/MRF/)

涉及的文件有：
- MRFZone.C|H
- MRFZones.C|H
- MRFZoneTemplates.C

其中of41的MRF文件更多一些，因为它把`inline function`内联函数给放到了`MRFZoneI.H`中了，并额外增加了`MRFZoneList.C|H`，`IOMRFZoneList.C|H`文件来处理`MRFZoneList`这个类。

### 代码编译

- 编译MRF模块

  - 需要修改OpenFOAM主体代码，并重新编译OpenFOAM的核心库`finiteVolume`。

  - 我觉得有点儿诡异，因为of41是自带MRF库的，所以打开`MRFZone.H`看看有什么特别之处。

  - ```c++
    // 查看of41, of2.0.x, AETK的MRFZone.H 之间的差异
    // #include文件
    //Common
    #include "dictionary.H"
    #include "wordList.H"
    #include "labelList.H"
    #include "dimensionedScalar.H"
    #include "dimensionedVector.H"
    #include "volFieldsFwd.H"
    #include "surfaceFields.H"
    #include "fvMatricesFwd.H"
    //AETKv1, of2.0.x
    #include "fvMatrices.H"
    //of41
    #include "mapPolyMesh.H"
    #include "Function1.H"
    #include "autoPtr.H"
    ```

  - 从`#include`文件来看，二者差不多，of41似乎还多一些功能。另外我注意到这个MRF库是Oliver Borm写的，而不是Chun Shen写的。因为Oliver Borm是DensityBasedTurbo的作者，所以可以初步推断AETK是基于DensityBasedTurbo开发的。而of41的MRF库不知道是谁写的。

  - 从Oliver Borm写的[源代码的说明](https://sourceforge.net/p/openfoam-extend/DensityBasedTurbo/ci/master/tree/src/)来看，它写MRF的目的是在OpenFOAM-1.6-ext中打补丁。所以很可能这个补丁已经不需要了。

- 因此可以直接去掉MRF编译

  - 然后发现EulerLocalDdtScheme有错误，考虑到of41中已经有了localEuler的ddt格式，先把它删除了先

  - ```shell
    cd src/timeStepping
    git rm EulerLocalDdtScheme -rf
    git rm MRF -rf
    vim Make/files # 去掉EulerLocalDdtScheme
    wmake
    ```

  - 提示错误：没有`dimensionedInternalField`这是由于of41的API变化引起的，具体可以参考[OpenFOAM®  Developer Upgrade Guide](http://www.openfoam.com/documentation/developer-upgrade-guide.php)，也比较容易更正：

  - ```shell
    rm lnInclude -rf
    find . -name "*.[C|H]" |xargs sed -i -e "s/dimensionedInternalField/internalField/g" 
    wmake
    ```

  - 此时又出现链接错误，找不到一些库。打开`Make/options`发现库的顺序和非常诡异，考虑到of41的API又有变化。所以这里需要大改。

  - 改成这样

  - ```shell
    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude  \
        -I$(LIB_SRC)/OpenFOAM/lnInclude  \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/specie/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/dynamicMesh/lnInclude 

    LIB_LIBS = \
        -lfiniteVolume \
        -lmeshTools \
        -lOpenFOAM \
        -ltriSurface \
        -lspecie \
        -lfiniteVolume \
        -lfvOptions \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lturbulenceModels \
        -lcompressibleTurbulenceModels 
    ```

  - 编译成功

- 编译godunovFlux模块

  - ```shell
    cd ../godunovFlux
    wmake
    ```

  - OK

- 编译inv模块，用于5*5矩阵求逆。

  - ```shell
    cd ../inv
    wmake
    ```

  - OK

- 编译求解器allSpeedUnsteadyFoam

  - 编译出错，错误信息如下：

  - ```shell
    Making dependency list for source file allSpeedUnsteadyFoam.C
    could not open file basicPsiThermo.H for source file allSpeedUnsteadyFoam.C due to No such file or directory
    could not open file turbulenceModel.H for source file allSpeedUnsteadyFoam.C due to No such file or directory
    could not open file MRFZones.H for source file allSpeedUnsteadyFoam.C due to No such file or directory
    g++ -std=c++0x -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I../godunovFlux -I../timeStepping/localTimeStep -I/opt/openfoam4/src/finiteVolume/lnInclude -I/opt/openfoam4/src/thermophysicalModels/specie/lnInclude -I/opt/openfoam4/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam4/src/turbulenceModels/compressible/turbulenceModel/lnInclude -I/opt/openfoam4/src/dynamicMesh/dynamicFvMesh/lnInclude -I/opt/openfoam4/src/dynamicMesh/dynamicMesh/lnInclude -I/opt/openfoam4/src/meshTools/lnInclude -I/opt/openfoam4/src/dynamicFvMesh/lnInclude  -IlnInclude -I. -I/opt/openfoam4/src/OpenFOAM/lnInclude -I/opt/openfoam4/src/OSspecific/POSIX/lnInclude   -fPIC -c allSpeedUnsteadyFoam.C -o Make/linux64GccDPInt32Opt/allSpeedUnsteadyFoam.o
    allSpeedUnsteadyFoam.C:59:28: fatal error: basicPsiThermo.H: No such file or directory
     #include "basicPsiThermo.H"
                                ^
    compilation terminated.
    /opt/openfoam4/wmake/rules/General/transform:8: recipe for target 'Make/linux64GccDPInt32Opt/allSpeedUnsteadyFoam.o' failed
    make: *** [Make/linux64GccDPInt32Opt/allSpeedUnsteadyFoam.o] Error 1
    ```

  - 从错误信息来看，找不到大量的文件，应该是`Make/options`有问题，参考rhoCentralDyMFoam的`Make/options`对其进行修改，如下：

  - ```shell
    EXE_INC = \
        -I../godunovFlux/lnInclude \
        -I../timeStepping/localTimeStep/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/dynamicFvMesh/lnInclude \
        -I$(LIB_SRC)/dynamicMesh/lnInclude \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/specie/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude

    EXE_LIBS = \
        -lGodunovFlux \
        -lDensityBasedSolvers \
        -lfiniteVolume \
        -ldynamicFvMesh \
        -ldynamicMesh \
        -ltopoChangerFvMesh \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lcompressibleTurbulenceModels \
        -lturbulenceModels \
        -lmeshTools
    ```

  - 编译还是出错，看来可能是api名称改变了。同样参考rhoCentralDyMFoam的源程序做修改。

    - `basicPsiThermo.H` =>`psiThermo.H` 
    - `turbulenceFluidThermoModel.H` =>`turbulentFluidThermoModel.H` 

  - 编译还是出错，显示找不到那堆SlopeLimiter。我就日了狗了，明明这堆东西就在那儿啊。

    - 仔细查看wmake的输出，发现画蛇添足地把lnInclude加到了godunovFlux后面，去掉就好了。
    - 但是其实这样不是很符合OpenFOAM规范的样子。

  - 编译终于出了正常的错误，也就是因为of41的API更新了。具体可以参考[OpenFOAM®  Developer Upgrade Guide](http://www.openfoam.com/documentation/developer-upgrade-guide.php)，也比较容易更正

    - 有的internalField应该是primitiveField，或者primitiveFieldRef

  - 编译出现连接错误，找不到`libGodunovFlux`和`libDensityBasedSolvers`，这肯定是`Make/options`文件的锅，在`EXE_LIBS`中加入` -L$(FOAM_USER_LIBBIN) `作为搜索路径。

### MRFZone详细对比
  - ```c++
    class MRFZone
    {
      //priviate data
      //Common
      const fvMesh& mesh_;
      const word name_;  
      const wordReList excludedPatchNames_;
      labelList excludedPatchLabels_;
      labelList internalFaces_;//- Internal faces that are part of MRF
      labelListList includedFaces_;//- Outside faces (per patch) that move with the MRF
      labelListList excludedFaces_;  //- Excluded faces (per patch) that do not move with the MRF
      
      //AETK
      const dictionary dict_;
      const dimensionedVector origin_;
      dimensionedVector axis_;
      const dimensionedScalar omega_;
      dimensionedVector Omega_;
      
      //of41
      dictionary coeffs_;
      const vector origin_;//- Origin of the axis
      vector axis_;  //- Axis vector
      autoPtr<Function1<scalar>> omega_;//- Angular velocty (rad/sec)
      bool active_;//- MRF region active flag
      word cellZoneName_;//- Name of cell zone
      label cellZoneID_;//- Cell zone ID
    }
    ```

  - 从代码来看，AETK版本的MRF私有数据比of41少，但是内容类型有些不同，功能应该差不多，of41额外多了`active`，`cellZone`数据。

  - ```c++
      class MRFZone
      {
        //private member function
        //common 
        //- Disallow default bitwise copy construct
        MRFZone(const MRFZone&);
        //- Disallow default bitwise assignment
        void operator=(const MRFZone&);
        void setMRFFaces();//- Divide faces in frame according to patch

        //AETKv1.0
        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType> void RelativeRhoFlux(const RhoFieldType& rho,surfaceScalarField& phi) const;
        //- Make the given relative mass/vol flux absolute within the MRF region
        template<class RhoFieldType> void AbsoluteRhoFlux(const RhoFieldType& rho,surfaceScalarField& phi) const;

        //of41
        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType> void makeRelativeRhoFlux(const RhoFieldType& rho,surfaceScalarField& phi) const;
        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType> void makeRelativeRhoFlux(const RhoFieldType& rho,FieldField<fvsPatchField, scalar>& phi) const;
        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType> void makeRelativeRhoFlux(const RhoFieldType& rho,Field<scalar>& phi,const label patchi) const;
        //- Make the given relative mass/vol flux absolute within the MRF region
        template<class RhoFieldType> void makeAbsoluteRhoFlux(const RhoFieldType& rho,surfaceScalarField& phi) const;
      }
      ```

  - 从私有函数来看，共有的是禁止默认的`copy construct`的部分和`setMRFFaces()`函数，其他方面of41多了几个函数，而且名字前加了`make`

  - ```c++
      class MRFZone
      {
        //common
        ClassName("MRFZone");
        MRFZone(const fvMesh& mesh, Istream& is);
        autoPtr<MRFZone> clone() const;
        void updateMesh(const mapPolyMesh& mpm);
        void correctBoundaryVelocity(volVectorField& U) const;
        
        //AETKv1.0
        class iNew{const fvMesh& mesh_;}
        void addCoriolis(fvVectorMatrix& UEqn) const;
        void addCoriolis(const volScalarField& rho, fvVectorMatrix& UEqn) const;
        void addCoriolis(const volScalarField& rho,const volVectorField& U,volVectorField& rhoUFlux) const;
        void relativeVelocity(const volVectorField& U,volVectorField& Urel) const;
        void relativeVelocity(volVectorField& U) const;
        void absoluteVelocity(volVectorField& U) const;
        void relativeFlux(surfaceScalarField& phi) const;
        void relativeFlux(const surfaceScalarField& rho,surfaceScalarField& phi) const;
        void absoluteFlux(surfaceScalarField& phi) const;
        void absoluteFlux(const surfaceScalarField& rho,surfaceScalarField& phi) const;
        void faceU(surfaceVectorField& zoneFaceU) const;//- Compute the pseudo face velocity of the MRF region
        void Su(const volScalarField& phi,const volVectorField& gradPhi,volScalarField& source) const;
        void Su(const volVectorField& phi,const volTensorField& gradPhi,volVectorField& source) const;
        void omega(volVectorField& source) const;
        friend Ostream& operator<<(Ostream& os, const MRFZone&);
        //of41
        void addCoriolis(fvVectorMatrix& UEqn,const bool rhs = false) const;
        void addCoriolis(const volScalarField& rho,fvVectorMatrix& UEqn,const bool rhs = false) const;
        void addCoriolis(const volVectorField& U,volVectorField& ddtU) const;
        void makeRelative(volVectorField& U) const;
        void makeRelative(surfaceScalarField& phi) const;
        void makeRelative(FieldField<fvsPatchField, scalar>& phi) const;
        void makeRelative(Field<scalar>& phi, const label patchi) const;
        void makeRelative(const surfaceScalarField& rho,surfaceScalarField& phi) const;
        void makeAbsolute(volVectorField& U) const;
        void makeAbsolute(surfaceScalarField& phi) const;
        void makeAbsolute(const surfaceScalarField& rho,surfaceScalarField& phi) const;
        void writeData(Ostream& os) const;
        bool read(const dictionary& dict);
      }
      ```

  - 从上述的MRF公开函数来看，MRF的功能都差不多，除了AETK多了`faceU`, `Su`, `omega()`函数之外，AETK和of41的MRF库主要是函数名字的差别。

  - 所以初步判断，of41和AETK的MRF库应该是可以兼容的，主要是api名称的变化而已；



# 参考文献

1. [DensityBasedTurbo的ppt](http://www.appliedccm.com/wp-content/uploads/2015/07/oliver_borm_slides.pdf)
2. [Implementation of density-based solver for all speeds in the framework of OpenFOAM](http://www.sciencedirect.com/science/article/pii/S0010465514002136)
