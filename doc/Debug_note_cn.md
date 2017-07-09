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

  - OK，allSpeedUnsteadyFoam编译成功

# 运行

AETK压缩包中自带了3个算例。

先从最简单的cavity开始，一运行就出错。

```
--> FOAM FATAL ERROR:
Unknown psiThermo type hPsiThermo<pureMixture<constTransport<specieThermo<hConstThermo<perfectGas>>>>>

Valid psiThermo types are:

44
(
hePsiThermo<homogeneousMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<homogeneousMixture<sutherland<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<homogeneousMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<inhomogeneousMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<inhomogeneousMixture<sutherland<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<inhomogeneousMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<multiComponentMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<multiComponentMixture<const<hConst<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<multiComponentMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<multiComponentMixture<sutherland<janaf<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<pureMixture<const<eConst<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<pureMixture<polynomial<hPolynomial<PengRobinsonGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<pureMixture<polynomial<janaf<PengRobinsonGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<pureMixture<sutherland<eConst<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<pureMixture<sutherland<hConst<PengRobinsonGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<pureMixture<sutherland<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<pureMixture<sutherland<hConst<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<pureMixture<sutherland<janaf<PengRobinsonGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<pureMixture<sutherland<janaf<PengRobinsonGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<pureMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<pureMixture<sutherland<janaf<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<reactingMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<reactingMixture<const<hConst<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<reactingMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<reactingMixture<sutherland<janaf<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<singleStepReactingMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<singleStepReactingMixture<sutherland<janaf<perfectGas<specie>>,sensibleInternalEnergy>>>
hePsiThermo<veryInhomogeneousMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<veryInhomogeneousMixture<sutherland<hConst<perfectGas<specie>>,sensibleEnthalpy>>>
hePsiThermo<veryInhomogeneousMixture<sutherland<janaf<perfectGas<specie>>,sensibleEnthalpy>>>
heheuPsiThermo<egrMixture<const<hConst<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<egrMixture<sutherland<janaf<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<egrMixture<sutherland<janaf<perfectGas<specie>>,absoluteInternalEnergy>>>
heheuPsiThermo<homogeneousMixture<const<hConst<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<homogeneousMixture<sutherland<janaf<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<homogeneousMixture<sutherland<janaf<perfectGas<specie>>,absoluteInternalEnergy>>>
heheuPsiThermo<inhomogeneousMixture<const<hConst<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<inhomogeneousMixture<sutherland<janaf<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<inhomogeneousMixture<sutherland<janaf<perfectGas<specie>>,absoluteInternalEnergy>>>
heheuPsiThermo<veryInhomogeneousMixture<const<hConst<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<veryInhomogeneousMixture<sutherland<janaf<perfectGas<specie>>,absoluteEnthalpy>>>
heheuPsiThermo<veryInhomogeneousMixture<sutherland<janaf<perfectGas<specie>>,absoluteInternalEnergy>>>
)

```

显然，要在上面的东西里找一个和` hPsiThermo<pureMixture<constTransport<specieThermo<hConstThermo<perfectGas>>>>>`相似的，应该是`hePsiThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>`，所以根据热物理包的生成关系，对`constant/thermophysicalProperties`进行修正。

```c++
//主要部分如下：
thermoType      hePsiThermo<pureMixture<const<hConst<perfectGas<specie>>,sensibleEnthalpy>>>;
mixture
{
    specie
    {
        nMoles          1;
        molWeight       28.9;
    }
    thermodynamics
    {
        Cp              1004.5;//1000;
        Hf              0;//2.544e+06;
    }
    transport
    {
        mu              1.72e-5;
        Pr              0.707;
    }
}
```

然后再次`foamJob allSpeedUnsteadyFoam`运行一次，出错：

```
--> FOAM FATAL IO ERROR:
keyword laplacian((thermo:mu+(0.6*(rho*nut))),k) is undefined in dictionary "/home/di/OpenFOAM/di-4.1/solvers/AETKv1/run/cavity/system/fvSchemes.laplacianSchemes"
```

这他妈的是哪儿来的？

首先确认一下不是AETK代码的锅

```shell
 $ find ../../src -name "*.[C|H]"|xargs grep "laplacian"
./allSpeedUnsteadyFoam/solveUnsteadyFluid.H:                    + fvc::laplacian(turbulence->alphaEff(), h)
./allSpeedUnsteadyFoam/solveUnsteadyFluid.H:                    + fvc::laplacian
./allSpeedUnsteadyFoam_dualtime/solveUnsteadyFluid.H:                    + fvc::laplacian(turbulence->alphaEff(), h)
./allSpeedUnsteadyFoam_dualtime/solveUnsteadyFluid.H:                    + fvc::laplacian
```

这时候我突然想起来，我还没有编译`allSpeedUnsteadyFoam_dualtime`，无所谓啦，开心就好。

改改改：

```c++
laplacianSchemes
{
  default         none;
  laplacian(alphaEff,h) Gauss linear corrected;
  //laplacian((mu+(0.6*mut)),k) Gauss linear corrected;
  laplacian((thermo:mu+(0.6*(rho*nut))),k) Gauss linear corrected;
  laplacian(thermo:alpha,h) Gauss linear corrected;
}
```

然后时间格式又出错了：

```
--> FOAM FATAL IO ERROR:
Unknown ddt scheme EulerLocal

Valid ddt schemes are :

8
(
CoEuler
CrankNicolson
Euler
SLTS
backward
bounded
localEuler
steadyState
)

```

对应的`fvSchemes`字段是：

```c++
ddtSchemes
{
    default none;

    ddt(rho)       EulerLocal  physDeltaT CoDeltaT;//backwardDual physDeltaT CoDeltaT;//
    ddt(rhoU)      EulerLocal  physDeltaT CoDeltaT;//backwardDual physDeltaT CoDeltaT;//
    ddt(rhoE)      EulerLocal  physDeltaT CoDeltaT;//backwardDual physDeltaT CoDeltaT;//
        ddt(p)             EulerLocal  physDeltaT CoDeltaT;//backwardDual physDeltaT CoDeltaT;//
        ddt(U)             EulerLocal  physDeltaT CoDeltaT;//backwardDual physDeltaT CoDeltaT;//
        ddt(T)             EulerLocal  physDeltaT CoDeltaT;//backwardDual physDeltaT CoDeltaT;//
        ddt(h)             EulerLocal  physDeltaT CoDeltaT;//backwardDual physDeltaT CoDeltaT;//
/*
    ddt(rho)       Euler;
    ddt(rhoU)      Euler;
    ddt(rhoE)      Euler;
        ddt(p)            Euler;
        ddt(U)            Euler;
        ddt(T)            Euler;
        ddt(h)            Euler;
*/
}
```

呵呵，我还有什么话说呢，真的猛士，敢于直面惨淡的人生，干预正视淋漓的鲜血。我不是猛士，我就用of41里有的试试吧，但是一看of有几种LTS，还是得斟酌一番：

- localEuler: 似乎是标准的LTS

  - LTS场从数据库（objectRegistry）中被读取出来。
  - LTS场的名字固定的，就是`rDeltaT`和`rSubDeltaT`。
  - 通常求解器中会包含`$FOAM_SRC/finiteVolume/cfdTools/general/include/createRDeltaT.H`来生成`rDeltaT`。
  - 从代码来看，这个`createRDeltaT.H`只能和`localEuler`配合使用。

- CoEuler: Courant数限制的LTS

  - 局部Courant数不超过指定值

  - 用法：

  - ```c++
    ddtSchemes
    {
      default SLTS phiName rhoName maxCo;
      // phiName: face flux
      // rhoName: density
      // maxCo: max local Courant number
    }
    ```

  - ​

- SLTS: 稳定化的LTS格式

  - 保证对流方程对角占优

  - 用法：

  - ```c++
    ddtSchemes
    {
      default SLTS phiName rhoName alpha;
      // phiName: face flux
      // rhoName: density
      // alpha: under relaxation factor
    }
    ```

然后试了以下这个选项：
```c++
ddtSchemes
{
  default localEuler;
}
```
出错了：
```
--> FOAM FATAL ERROR:

    request for volScalarField rSubDeltaT from objectRegistry region0 failed
    available objects of type volScalarField are

26
(
thermo:mu
-div(rhoFlux)
(-div(rhoUFlux)+div(tau)).component(2)
thermo:psi
U.component(0)
(1|(((((Cp|Cv)*(Cp-Cv))*T)*(((1|Ur2)-(1|(((Cp|Cv)*(Cp-Cv))*T)))+(1|((Cp-Cv)*T))))-((Cp|Cv)-1)))
h
U.component(1)
rho
(((((Cp|Cv)-1)*(1|(((((Cp|Cv)*(Cp-Cv))*T)*(((1|Ur2)-(1|(((Cp|Cv)*(Cp-Cv))*T)))+(1|((Cp-Cv)*T))))-((Cp|Cv)-1))))*T)|rho)
TOL
(-div(rhoUFlux)+div(tau)).component(1)
magSqr(U)
thermo:psi_0
rhoEFlux2
p
T
rhoE
(((-div(rhoEFlux)+div((tau&U)))+laplacian(interpolate(thermo:alpha),h))+laplacian(interpolate((thermo:mu+(0.6*(rho*nut)))),k))
(((Cp|Cv)-1)*(1|(((((Cp|Cv)*(Cp-Cv))*T)*(((1|Ur2)-(1|(((Cp|Cv)*(Cp-Cv))*T)))+(1|((Cp-Cv)*T))))-((Cp|Cv)-1))))
Ur2
U.component(2)
(-div(rhoUFlux)+div(tau)).component(0)
(((1|Ur2)-(1|(((Cp|Cv)*(Cp-Cv))*T)))+(1|((Cp-Cv)*T)))
rhoFlux2
thermo:alpha
)


    From function const Type& Foam::objectRegistry::lookupObject(const Foam::word&) const [with Type = Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh>]
    in file /opt/openfoam4/src/OpenFOAM/lnInclude/objectRegistryTemplates.C at line 193.

FOAM aborting

#0  Foam::error::printStack(Foam::Ostream&) at ??:?
#1  Foam::error::abort() at ??:?
#2  Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> const& Foam::objectRegistry::lookupObject<Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> >(Foam::word const&) const at ??:?
#3  Foam::fv::localEulerDdtScheme<double>::fvmDdt(Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> const&) at ??:?
#4  Foam::tmp<Foam::fvMatrix<double> > Foam::fvm::ddt<double>(Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> const&) at ??:?
#5  ? at ??:?
#6  __libc_start_main in "/lib/x86_64-linux-gnu/libc.so.6"
#7  ? at ??:?
```

艹，什么乱七八糟的。

仔细看EulerLocal的代码，发现这个程序自己实现的EulerLocal同时包含物理时间步和LTS时间步，在fvSchemes中需要指定它俩的名称。即算例中的`PhysDeltaT, CoDeltaT`，不明觉厉，似乎用localEuler代替EulerLocal有些障碍。

那就把EulerLocal编译一下试试。将`EulerLocalDdtScheme.C`加入`Make/files`

编译出错：

```
EulerLocalDdtScheme/EulerLocalDdtScheme.C:52:1: error: redefinition of ‘Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> > Foam::fv::EulerLocalDdtScheme<Type>::fvcDdt(const Foam::dimensioned<Type>&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:118:1: error: redefinition of ‘Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> > Foam::fv::EulerLocalDdtScheme<Type>::fvcDdt(const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:176:1: error: redefinition of ‘Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> > Foam::fv::EulerLocalDdtScheme<Type>::fvcDdt(const dimensionedScalar&, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:235:1: error: redefinition of ‘Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> > Foam::fv::EulerLocalDdtScheme<Type>::fvcDdt(const volScalarField&, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:297:1: error: redefinition of ‘Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::EulerLocalDdtScheme<Type>::fvmDdt(const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:345:1: error: redefinition of ‘Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::EulerLocalDdtScheme<Type>::fvmDdt(const dimensionedScalar&, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:389:1: error: redefinition of ‘Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::EulerLocalDdtScheme<Type>::fvmDdt(const volScalarField&, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:435:1: error: redefinition of ‘Foam::tmp<typename Foam::fv::ddtScheme<Type>::fluxFieldType> Foam::fv::EulerLocalDdtScheme<Type>::fvcDdtPhiCorr(const volScalarField&, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&, const fluxFieldType&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:495:1: error: redefinition of ‘Foam::tmp<typename Foam::fv::ddtScheme<Type>::fluxFieldType> Foam::fv::EulerLocalDdtScheme<Type>::fvcDdtPhiCorr(const volScalarField&, const volScalarField&, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&, const fluxFieldType&)’
EulerLocalDdtScheme/EulerLocalDdtScheme.C:624:25: error: redefinition of ‘Foam::tmp<Foam::GeometricField<double, Foam::fvsPatchField, Foam::surfaceMesh> > Foam::fv::EulerLocalDdtScheme<Type>::meshPhi(const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&)’
```

全是`redifinition`错误，这是个啥啊，具体点进去看各个文件出错的地方，大概是这样的：

- 一般类分为A.H和A.C文件，A.C会包含A.H，但是A.H不会包含A.C了。

- 奇葩的EulerLocal（包括OpenFOAM）在A.H的末尾加上了NoRepositiory的宏定义，如果有这个宏（通过`gcc -DNoRepositiory`传入），那么A.H会把A.C包含进来。

  - 关于NoRepository的问题，在http://xiaopingqiu.github.io/2016/04/25/TurbulenceModel-30-macro/ 有所讨论。
  - 和类模板的实例化有关。

- 然后A.C包含A.H包含A.C，这个事儿就乱了套了。

- 这种重复定义有两种解决方法，参考：https://stackoverflow.com/a/707954/4592964

  - `#pragma once`，非标准，但是被vc和gcc支持。

  - ````c++
    #ifndef XXX
    #define XXX
    //...
    #endif
    ````

  - 但是，似乎加了还是没什么效果。

- 仔细查看localEuler的定义，发现在`Make/files`中包含的不是`localEulerDdtScheme.H`，而是`localEulerDdtSchemes.H`，多的一个s表示static静态的意思，好象是要将其实例化。

  - 果然，重复定义的问题解决了。但是又出来了其他的问题：

- ```
  In file included from /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionSet.H:45:0,
                   from /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionedType.H:40,
                   from /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:39,
                   from EulerLocalDdtScheme/EulerLocalDdtScheme.H:50,
                   from EulerLocalDdtScheme/EulerLocalDdtSchemes.C:30:
  /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:260:53: error: expected initializer before ‘<’ token
       defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                       ^
  /opt/openfoam4/src/OpenFOAM/lnInclude/className.H:91:24: note: in definition of macro ‘defineTypeNameWithName’
       const ::Foam::word Type::typeName(Name)
                          ^
  /opt/openfoam4/src/OpenFOAM/lnInclude/className.H:112:5: note: in expansion of macro ‘defineTemplateTypeNameWithName’
       defineTemplateTypeNameWithName(Type, Type::typeName_())
       ^
  /opt/openfoam4/src/OpenFOAM/lnInclude/className.H:135:5: note: in expansion of macro ‘defineNamedTemplateTypeName’
       defineNamedTemplateTypeName(Type);                                         \
       ^
  /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:260:5: note: in expansion of macro ‘defineNamedTemplateTypeNameAndDebug’
       defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
       ^
  /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:274:1: note: in expansion of macro ‘makeFvDdtTypeScheme’
   makeFvDdtTypeScheme(SS, vector)                                                \
   ^
  EulerLocalDdtScheme/EulerLocalDdtSchemes.C:39:5: note: in expansion of macro ‘makeFvDdtScheme’
       makeFvDdtScheme(EulerLocalDdtScheme)^M
       ^
  ```

- 说一下gcc的错误信息：

  - ```
    In file included from /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionSet.H:45:0,
                     from /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionedType.H:40,
                     from /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:39,
                     from EulerLocalDdtScheme/EulerLocalDdtScheme.H:50,
                     from EulerLocalDdtScheme/EulerLocalDdtSchemes.C:30:
    ```

  - 上面说的是文件的include关系，下面的在指定的行include 上面的，第一行不知是啥，但是最后一行是正在编译的文件。验证如下：

  - ```shell
    $ sed -n "6,10p" log
    In file included from /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionSet.H:45:0,
                     from /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionedType.H:40,
                     from /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:39,
                     from EulerLocalDdtScheme/EulerLocalDdtScheme.H:50,
                     from EulerLocalDdtScheme/EulerLocalDdtSchemes.C:30:
    $ sed -n "30p" EulerLocalDdtScheme/EulerLocalDdtSchemes.C
    #include "EulerLocalDdtScheme.H"
    $ sed -n "50p" EulerLocalDdtScheme/EulerLocalDdtScheme.H
    #include "ddtScheme.H"
    $ sed -n "39p" /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H
    #include "dimensionedType.H"
    $ sed -n "40p" /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionedType.H
    #include "dimensionSet.H"
    $ sed -n 45p /opt/openfoam4/src/OpenFOAM/lnInclude/dimensionSet.H
    #include "className.H"
    ```

  - 后面的信息应该是错误信息`error`和调用关系`note`

  - ```
    /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:260:53: error: expected initializer before ‘<’ token
         defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                         ^
    /opt/openfoam4/src/OpenFOAM/lnInclude/className.H:91:24: note: in definition of macro ‘defineTypeNameWithName’
         const ::Foam::word Type::typeName(Name)
                            ^
    /opt/openfoam4/src/OpenFOAM/lnInclude/className.H:112:5: note: in expansion of macro ‘defineTemplateTypeNameWithName’
         defineTemplateTypeNameWithName(Type, Type::typeName_())
         ^
    /opt/openfoam4/src/OpenFOAM/lnInclude/className.H:135:5: note: in expansion of macro ‘defineNamedTemplateTypeName’
         defineNamedTemplateTypeName(Type);                                         \
         ^
    /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:260:5: note: in expansion of macro ‘defineNamedTemplateTypeNameAndDebug’
         defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
         ^
    /opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H:274:1: note: in expansion of macro ‘makeFvDdtTypeScheme’
     makeFvDdtTypeScheme(SS, vector)                                                \
     ^
    EulerLocalDdtScheme/EulerLocalDdtSchemes.C:39:5: note: in expansion of macro ‘makeFvDdtScheme’
         makeFvDdtScheme(EulerLocalDdtScheme)^M
         ^
    ```

  - 上面的也是从下网上溯源的：

    - `EulerLocalDdtSchemes.C`在第39行调用了宏`makeFvDdtScheme`

      ​

    - 宏`makeFvDdtScheme`是在`/opt/openfoam4/src/finiteVolume/lnInclude/ddtScheme.H`的第274行定义的

    - ```c++
      #define makeFvDdtScheme(SS)                                                    \
                                                                                     \
      makeFvDdtTypeScheme(SS, scalar)                                                \
      makeFvDdtTypeScheme(SS, vector)                                                \
      makeFvDdtTypeScheme(SS, sphericalTensor)                                       \
      makeFvDdtTypeScheme(SS, symmTensor)                                            \
      makeFvDdtTypeScheme(SS, tensor)                                                \
                                                                                     \
      namespace Foam                                                                 \
      {                                                                              \
      namespace fv                                                                   \
      {                                                                              \
                                                                                     \
      template<>                                                                     \
      tmp<surfaceScalarField> SS<scalar>::fvcDdtUfCorr                               \
      (                                                                              \
          const volScalarField& U,                                                   \
          const surfaceScalarField& Uf                                               \
      )                                                                              \
      {                                                                              \
          NotImplemented;                                                            \
          return surfaceScalarField::null();                                         \
      }                                                                              \
                                                                                     \
      template<>                                                                     \
      tmp<surfaceScalarField> SS<scalar>::fvcDdtPhiCorr                              \
      (                                                                              \
          const volScalarField& U,                                                   \
          const surfaceScalarField& phi                                              \
      )                                                                              \
      {                                                                              \
          NotImplemented;                                                            \
          return surfaceScalarField::null();                                         \
      }                                                                              \
                                                                                     \
      template<>                                                                     \
      tmp<surfaceScalarField> SS<scalar>::fvcDdtUfCorr                               \
      (                                                                              \
          const volScalarField& rho,                                                 \
          const volScalarField& U,                                                   \
          const surfaceScalarField& Uf                                               \
      )                                                                              \
      {                                                                              \
          NotImplemented;                                                            \
          return surfaceScalarField::null();                                         \
      }                                                                              \
                                                                                     \
      template<>                                                                     \
      tmp<surfaceScalarField> SS<scalar>::fvcDdtPhiCorr                              \
      (                                                                              \
          const volScalarField& rho,                                                 \
          const volScalarField& U,                                                   \
          const surfaceScalarField& phi                                              \
      )                                                                              \
      {                                                                              \
          NotImplemented;                                                            \
          return surfaceScalarField::null();                                         \
      }                                                                              \
                                                                                     \
      }                                                                              \
      }
      ```

    - 其中`SS`是宏的参数，行末的`\`字符是续行符；

    - 主要功能是展开成5个不同变量种类（Type=scalar/vector/tensor/sphericalTensor/symmTensor）的ddtScheme，然后把用scalar特化相关的一些无意义函数定义为无实现`NotImplemented`。

    - 然后追溯到两个参数的宏`makeFvDdtTypeScheme(SS, Type)`的定义，这个宏名字和单个参数的一样

    - ```c++
      #define makeFvDdtTypeScheme(SS, Type)                                          \
          defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                                     \
          namespace Foam                                                             \
          {                                                                          \
              namespace fv                                                           \
              {                                                                      \
                  ddtScheme<Type>::addIstreamConstructorToTable<SS<Type>>            \
                      add##SS##Type##IstreamConstructorToTable_;                     \
              }                                                                      \
          }
      ```

    - 这里调用了` defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0)`，特化为`defineNamedTemplateTypeNameAndDebug(Foam::fv::EulerLocalDdtScheme<Foam::vector>,0)`再次追溯：

    - ```c++
      //$vim ./OpenFOAM/db/typeInfo/className.H +135
      //- Define the typeName and debug information for templates
      #define defineNamedTemplateTypeNameAndDebug(Type, DebugSwitch)                 \
          defineNamedTemplateTypeName(Type);                                         \
          defineNamedTemplateDebugSwitch(Type, DebugSwitch)
      ```

    - 可以看出，0只是DebugSwitch而已，但是`Type=Foam::fv::EulerLocalDdtScheme<Foam::vector>`继续追溯：

    - ```c++
      //- Define the typeName directly for template classes
      #define defineNamedTemplateTypeName(Type)                                      \
          defineTemplateTypeNameWithName(Type, Type::typeName_())
      ```

    - 这显然又是一个宏的调用，继续追溯：

    - ```c++
      //- Define the typeName as \a Name for template classes
      #define defineTemplateTypeNameWithName(Type, Name)                             \
          template<>                                                                 \
          defineTypeNameWithName(Type, Name)
      ```

    - 继续追溯：

    - ```c++
      //- Define the typeName, with alternative lookup as \a Name
      #define defineTypeNameWithName(Type, Name)                                     \
          const ::Foam::word Type::typeName(Name)
      ```

    - 到了这儿，我TMD还是没看懂。直接在`Make/options`中加入`-E`选项，然后找出实际用的编译命令：

    - ```shell
      $ g++ -std=c++0x -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -E -I/opt/openfoam4/src/finiteVolume/lnInclude -I/opt/openfoam4/src/meshTools/lnInclude -I/opt/openfoam4/src/OpenFOAM/lnInclude -I/opt/openfoam4/src/transportModels/compressible/lnInclude -I/opt/openfoam4/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam4/src/thermophysicalModels/specie/lnInclude -I/opt/openfoam4/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam4/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam4/src/dynamicMesh/lnInclude  -IlnInclude -I. -I/opt/openfoam4/src/OpenFOAM/lnInclude -I/opt/openfoam4/src/OSspecific/POSIX/lnInclude   -fPIC -c EulerLocalDdtScheme/EulerLocalDdtSchemes.C -o file_cpp
      ```

    - 其中在`file_cpp`中可以找到上述代码在处理之后的样子是：

    - ```c++
      template<> const ::Foam::word Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName(Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName_()); template<> int Foam::fv::EulerLocalDdtScheme<Foam::scalar>::debug(::Foam::debug::debugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName_(), 0)); template<> const Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::scalar> > Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::scalar> >::registerDebugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName_()); namespace Foam { namespace fv { ddtScheme<scalar>::addIstreamConstructorToTable<EulerLocalDdtScheme<scalar>> addEulerLocalDdtSchemescalarIstreamConstructorToTable_; } } template<> const ::Foam::word Foam::fv::EulerLocalDdtScheme<Foam::vector>::typeName(Foam::fv::EulerLocalDdtScheme<Foam::vector>::typeName_()); template<> int Foam::fv::EulerLocalDdtScheme<Foam::vector>::debug(::Foam::debug::debugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::vector>::typeName_(), 0)); template<> const Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::vector> > Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::vector> >::registerDebugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::vector>::typeName_()); namespace Foam { namespace fv { ddtScheme<vector>::addIstreamConstructorToTable<EulerLocalDdtScheme<vector>> addEulerLocalDdtSchemevectorIstreamConstructorToTable_; } } template<> const ::Foam::word Foam::fv::EulerLocalDdtScheme<Foam::sphericalTensor>::typeName(Foam::fv::EulerLocalDdtScheme<Foam::sphericalTensor>::typeName_()); template<> int Foam::fv::EulerLocalDdtScheme<Foam::sphericalTensor>::debug(::Foam::debug::debugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::sphericalTensor>::typeName_(), 0)); template<> const Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::sphericalTensor> > Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::sphericalTensor> >::registerDebugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::sphericalTensor>::typeName_()); namespace Foam { namespace fv { ddtScheme<sphericalTensor>::addIstreamConstructorToTable<EulerLocalDdtScheme<sphericalTensor>> addEulerLocalDdtSchemesphericalTensorIstreamConstructorToTable_; } } template<> const ::Foam::word Foam::fv::EulerLocalDdtScheme<Foam::symmTensor>::typeName(Foam::fv::EulerLocalDdtScheme<Foam::symmTensor>::typeName_()); template<> int Foam::fv::EulerLocalDdtScheme<Foam::symmTensor>::debug(::Foam::debug::debugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::symmTensor>::typeName_(), 0)); template<> const Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::symmTensor> > Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::symmTensor> >::registerDebugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::symmTensor>::typeName_()); namespace Foam { namespace fv { ddtScheme<symmTensor>::addIstreamConstructorToTable<EulerLocalDdtScheme<symmTensor>> addEulerLocalDdtSchemesymmTensorIstreamConstructorToTable_; } } template<> const ::Foam::word Foam::fv::EulerLocalDdtScheme<Foam::tensor>::typeName(Foam::fv::EulerLocalDdtScheme<Foam::tensor>::typeName_()); template<> int Foam::fv::EulerLocalDdtScheme<Foam::tensor>::debug(::Foam::debug::debugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::tensor>::typeName_(), 0)); template<> const Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::tensor> > Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::tensor> >::registerDebugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::tensor>::typeName_()); namespace Foam { namespace fv { ddtScheme<tensor>::addIstreamConstructorToTable<EulerLocalDdtScheme<tensor>> addEulerLocalDdtSchemetensorIstreamConstructorToTable_; } } namespace Foam { namespace fv { template<> tmp<surfaceScalarField> EulerLocalDdtScheme<scalar>::fvcDdtUfCorr ( const volScalarField& U, const surfaceScalarField& Uf ) { ::Foam::FatalError((__PRETTY_FUNCTION__), "EulerLocalDdtScheme/EulerLocalDdtSchemes.C", 39) << "Not implemented" << ::Foam::abort(FatalError);; return surfaceScalarField::null(); } template<> tmp<surfaceScalarField> EulerLocalDdtScheme<scalar>::fvcDdtPhiCorr ( const volScalarField& U, const surfaceScalarField& phi ) { ::Foam::FatalError((__PRETTY_FUNCTION__), "EulerLocalDdtScheme/EulerLocalDdtSchemes.C", 39) << "Not implemented" << ::Foam::abort(FatalError);; return surfaceScalarField::null(); } template<> tmp<surfaceScalarField> EulerLocalDdtScheme<scalar>::fvcDdtUfCorr ( const volScalarField& rho, const volScalarField& U, const surfaceScalarField& Uf ) { ::Foam::FatalError((__PRETTY_FUNCTION__), "EulerLocalDdtScheme/EulerLocalDdtSchemes.C", 39) << "Not implemented" << ::Foam::abort(FatalError);; return surfaceScalarField::null(); } template<> tmp<surfaceScalarField> EulerLocalDdtScheme<scalar>::fvcDdtPhiCorr ( const volScalarField& rho, const volScalarField& U, const surfaceScalarField& phi ) { ::Foam::FatalError((__PRETTY_FUNCTION__), "EulerLocalDdtScheme/EulerLocalDdtSchemes.C", 39) << "Not implemented" << ::Foam::abort(FatalError);; return surfaceScalarField::null(); } } }
      ```

    - 提取出数据是scalar那一段是：

    - ```c++
      template<> const ::Foam::word Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName(Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName_());
      template<> int Foam::fv::EulerLocalDdtScheme<Foam::scalar>::debug(::Foam::debug::debugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName_(), 0));
      template<> const Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::scalar> > Foam::RegisterDebugSwitch<Foam::fv::EulerLocalDdtScheme<Foam::scalar> >::registerDebugSwitch(Foam::fv::EulerLocalDdtScheme<Foam::scalar>::typeName_());
      namespace Foam
      {
        namespace fv
        {
          ddtScheme<scalar>::addIstreamConstructorToTable<EulerLocalDdtScheme<scalar>> addEulerLocalDdtSchemescalarIstreamConstructorToTable_;
        }
      }
      ```

    - 还是没闹特别明白，经过搜索这个错误，我发现它的错误意思是在`<`这个符号(`token`)位置之前，需要有`EulerLocalDdtScheme`的初始化函数（定义？实现？）

    - 根据https://stackoverflow.com/questions/21997706/error-expected-initializer-before-using-c这个帖子，似乎有可能是某些地方缺少了`;`导致的！查一下代码先。

    - 对比了一下localEuler和EulerLocal的头文件，发现EulerLocal少了一个(fvMesh)的构造函数，加上试试先。发现无效。

    - 继续查找stack overflow，发现这个帖子提到了`does not name a type`可能是`namespace`相关的问题。但是没能具体定位。

    - 继续对比localEuler和EulerLocal的头文件，发现EulerLocal少了两个`fvcDdtUfCorr`函数;

      - 这有可能是原因，因为这个函数是在父虚类ddtScheme中定义为虚函数的，如果他们没有被实现，则EulerLocal还是一个虚基类，没有办法被实例化。可能就会导致前面的找不到initiator的问题。

      - perfect story，完美解释了所有现象。

      - 在EulerLocalDdtScheme.H中加入：

      - ```c++
        //only display codes to add
        class EulerLocalDdtScheme
        {
          public:
          ...
          tmp<fluxFieldType> fvcDdtUfCorr
            (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
          );
          tmp<fluxFieldType> fvcDdtUfCorr
            (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
          );
          ...
        }

        ...
        template<>
        tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtUfCorr
        (
            const GeometricField<scalar, fvPatchField, volMesh>& U,
            const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
        );
        template<>
        tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtUfCorr
        (
            const volScalarField& rho,
            const volScalarField& U,
            const surfaceScalarField& Uf
        );
        ```

      - 话说这个`fvcDdtUfCorr`函数还是of40引入的，似乎是为了配合`dotInterpolate(surfaceF, volF.oldTime)`函数（参考https://github.com/OpenFOAM/OpenFOAM-dev/commit/d52546b3ce776dbfcfc24eb3ffc30e80064d3b4c），返回的是`dotInterpolate`和`(surfaceF & surfaceF.oldTime)`的差相关的一个量。

      - 参考不同版本之间`pisoFoam`中`pEqn`中的修正项的差异可以发现：

        - of2x版本中是`fvc::ddtPhiCorr(rUA, U, phi)`，有三个参数（可压缩流时会加密度rho，形成4个参数）。
        - of2x中的`ddtPhiCorr(rUA,U,phi)`对应的项在of41的pEqn中变成了`fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)`
        - 而of41中的`fvc::ddtCorr`，仔细查看代码`$FOAM_SRC/finiteVolume/finiteVolume/fvc/fvcDdt.C +226`可以发现它们调用了`ddtScheme`中的`fvcDdtUfCorr()`函数。

      - 这涉及到面流量场和体速度场的相容性的问题，即插值后的U和phi是否一致。由于在很多分离算法（SIMPLE, PISO等）中，再考虑在同位网格条件下，phi和U由于更新有先后顺序之分，可能会导致U和phi不一致，所以需要修正。

      - 参考资料：https://zh.scribd.com/doc/48195039/ddtPhiCorr，对于of2x版本来说，以简单的不可压流为例：

        - `fvc::ddtPhiCorr(...)`调用`ddtScheme::fvcDdtPhiCorr(rA,U,phi)`，其计算公式如下：

        - $$
          \delta \phi=\phi^{o}-U^o\cdot S_f\\
          K_c = 1-\min{(1,\frac{|\delta\phi|}{|\phi^o|+\epsilon})}\\
          \phi_{c}=\delta \phi\cdot K_c\cdot \frac 1 A\cdot \frac 1 {\Delta t}
          $$

        - `ddtScheme::DdtPhiCorrCoeff()`用于计算$K_c$

        - 返回$\phi_c$

      - 对于of41而言，情况有所变化：

        - 没有了`fvc::ddtPhiCorr()`，有`fvc::ddtCorr()`，不再包含`rA`这一项参数。
        - `fvc::ddtCorr()`调用的函数是`ddtScheme::fvcDdtPhiCorr(U,phi)`或者`ddtScheme::fvcDdtUfCorr(U,Uf)`，根据调用参数的选择不同。
        - 因为`Uf`必须是矢量，`phi`是标量，所以后面会限制一下模板的特化。
        - 由于没有再包含`rA`这项参数，所以在很多算法中的`pEqn.H`中需要将of2x中的`fvc::ddtPhiCorr(rUA, U, phi)`更换为of41中的`fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)`

      - 所以EulerLocal的移植还真挺难的。真想重写。

      - 重写得看懂EulerLocal的语义。然后在of41的localEuler或者CoEuler的基础上进行修改。

      - 还是先试试改吧：`git checkout -b "fuckEulerLocal"`开个分支

      - 然后把`EulerLocalDdtScheme`类的`fvcDdtPhiCorr`函数的`rA`参数都去掉，相应的修改很琐碎。

      - 中途搜索了一下AETK的代码，没有用到`fvc::ddtPhiCorr`这个函数，那其他代码应该不用修改了。

      - 参考着localEuler的实现，把`fvcDdtUfCorr()`函数的实现给加到`EulerLocalDdtScheme.C`中去。

      - 重要的修改点：

      - ```c++
        //const surfaceScalarField rDeltaT(fvc::interpolate(localRDeltaT()));
        //上面替换成下面
        const objectRegistry& registry = this->mesh();
        // get access to the scalar beta[i]
        const scalarField& beta =
          registry.lookupObject<scalarField>(deltaTName_);
        volScalarField rDeltaT =
          1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));
        ```

      - 然后悲剧地发现，还是有` error: invalid new-expression of abstract class type ‘Foam::fv::EulerLocalDdtScheme<>`这样的错误。

      - 对比一下，发现of41的api改变太大了，`ddtScheme`的子类需要实现比of2x版本更多纯虚函数。包括：

      - | 纯虚函数个数        | of41 | of2x |
        | ------------- | ---- | ---- |
        | fvcDdt        | 5    | 4    |
        | fvmDdt        | 4    | 3    |
        | fvcDdtPhiCorr | 2    | 2    |
        | fvcDdtUfCorr  | 2    | 0    |
        | meshPhi       | 1    | 1    |
        | type          | 1    | 1    |

      - 还有以下函数需要实现

      - ```c++
        fvcDdt(const volScalarField &alpha, const volScalarField &rho, const GeometricField< Type, fvPatchField, volMesh > &)=0;
        fvmDdt(const volScalarField &alpha, const volScalarField &rho, const GeometricField< Type, fvPatchField, volMesh > &vf)=0;
        ```

      - 偷个懒，尝试用`NotImplemented`宏试试：

      - ```c++
        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
          (
          const volScalarField& alpha,
          const volScalarField& rho,
          const GeometricField<Type, fvPatchField, volMesh>& psi
        );
        template<class Type>
          tmp<GeometricField<Type, fvPatchField, volMesh>>
            localEulerDdtScheme<Type>::fvcDdt
              (
              const volScalarField& alpha,
              const volScalarField& rho,
              const GeometricField<Type, fvPatchField, volMesh>& vf
            )
            {
              NotImplemented;
              return tmp<GeometricField<Type, fvPatchField, volMesh>>();//return null tmp
            }

        tmp<fvMatrix<Type>> fvmDdt
          (
          const volScalarField& alpha,
          const volScalarField& rho,
          const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        template<class Type>
          tmp<fvMatrix<Type>>
          EulerLocalDdtScheme<Type>::fvmDdt
            (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& vf
          )
          {
            NotImplemented;
            return tmp<fvMatrix<Type>>();//return null tmp
          }
        ```

      - Bingo!终于出了正常错误。61个。。。但是主要的错误类型是：

      - ```c++
        //error: fvcDdtPhiCoeff not in the scope
        //cause: api change, or c++ change ...
        this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime());
        //solution:
        //:s/fvcDdtPhiCoeff/this->fvcDdtPhiCoeff/g
        //result: error from 61->45

        //error: binding ‘A’ to reference of type ‘B’ discards qualifiers
        //cause: tmp<> API change;
        fvMatrix<Type>& fvm = tfvm();
        //solution: 
        //grep "tfvm()" EulerLocalDdtScheme/EulerLocalDdtScheme.C
        //sed -i -e "s/tfvm()/tfvm.ref()/g" EulerLocalDdtScheme/EulerLocalDdtScheme.C
        //result: error from 45->30

        //error: expected initializer before ‘<’ token
        //cause: makeFvDdtScheme(EulerLocalDdtScheme) is enclosed by Foam::fv namespace
        //solution: remove the Foam::fv namespace in EulerLocalDdtSchemes.C
        //result: No. of error reduced to 12

        //error: fvcDdtPhiCoeff no matching problem
        //cause: not invoke fvcDdtPhiCoeff(vF,sF[,sF2]) correctly.
        //diagnostics: using tricks to find the type of each part of the arguments from error info from compiler
        //int a=0;
        //a = fvc::interpolate(rDeltaT)*phi.oldTime();
        //a = fvc::interpolate(rDeltaT*U.oldTime())&mesh().Sf();
        //a = fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phi.oldTime());
        //output shows fvcDdtPhiCoeff is source of error; looking into it I found the `rho` should not be there 
        //independently.
        //solution: invoke it this way: fvcDdtPhiCoeff(rho.oldTime()*U.oldTime(),phi.oldTime())
        //result: reduction to 8 errors.

        //error: no matching function for call to ‘Foam::GeometricField<Foam::Vector<double>, 
        //Foam::fvsPatchField, Foam::surfaceMesh>::GeometricField(Foam::IOobject, 
        //Foam::tmp<Foam::Field<Foam::Vector<double> > >)’
        //cause: vF*sF result in common superclass type Field rather than either vF or sF
        //solutioin: 
        // change rDeltaT to fvc::interpolate(rDeltaT)
        // result: 0 errors.
        ```

      - 终于成功了，因为是动态链接库，所以不用重新生成程序。直接就可以用。

  - 如果需要动态调试模板，可以用gcc 的 `-gdwarf-3/4`选项配合gdb的`marco`命令。

  尝试了一下，cavity是可以运行了，但是似乎只有`Ux, Uy`的信息，会不会只能计算二维算例呢？看了一眼代码，不是这样的。

  但是运行并不算很快的样子。

  尝试并行：OK！

  打开log查看，发现有`Duplicate entry`的错误。

  ```
  Duplicate entry limitedVanLeer in runtime selection table limitedSurfaceInterpolationScheme
  #0      /opt/openfoam4/platforms/linux64GccDPInt32Opt/lib/libOpenFOAM.so(_ZN4Foam5error14safePrintStackERSo+0x3e) [0x7fe9119042ee]
  #1      /home/di/OpenFOAM/di-4.1/platforms/linux64GccDPInt32Opt/lib/libGodunovFlux.so(+0x76d59) [0x7fe916666d59]
  #2      /lib64/ld-linux-x86-64.so.2(+0x104ea) [0x7fe9169ba4ea]
  #3      /lib64/ld-linux-x86-64.so.2(+0x105fb) [0x7fe9169ba5fb]
  #4      /lib64/ld-linux-x86-64.so.2(+0xcfa) [0x7fe9169aacfa]
  Duplicate entry limitedVanLeer in runtime selection table limitedSurfaceInterpolationScheme
  Duplicate entry limitedVanLeer in runtime selection table surfaceInterpolationScheme
  ...
  ```

  果断用官方的呀。删除vanLeer这个slopeLimiter。

## 编译allSpeedUnsteadyFoam_dualTime

除了相同的问题外，这个求解器还有一些独特的问题：

```c++
// error: ‘velocity_Inlet’ was not declared in this scope
// cause: 作者代码貌似不全，这个velocity_Inlet可全局搜索，是定义在system/controlDict中的，意义是文献中全局固定的参考速度。
// solution: loop前加入
scalar velocity_Inlet = runTime.controlDict().lookupOrDefault<scalar>("velocity_Inlet", 3.0);

// error: passing ‘const Internal {aka const Foam::DimensionedField<double, Foam::volMesh>}’ as ‘this’ argument discards qualifiers [-fpermissive]
// error: no matching function for call to ‘gMax(Foam::tmp<Foam::DimensionedField<double, Foam::volMesh> >)’
// error: no match for ‘operator=’ (operand types are ‘const Internal {aka const Foam::DimensionedField<double, Foam::volMesh>}’ and ‘Foam::tmp<Foam::Field<double> >’)
// cause: api 升级
// solution: 替换internalField()为primitiveField()或 primitiveFieldRef()
```

然后编译是可以了，但是在运行`forwardStep`这个算例的时候遇到障碍。

```
--> FOAM FATAL ERROR:

    request for volScalarField rhoOld from objectRegistry region0 failed
    available objects of type volScalarField are

21
(
rhoEFlux3
tmp_velocity
thermo:mu
tmp_zero
rhoOldOld
thermo:psi
rhoE_0
h
rho
rhoEOldOld
TOL
thermo:psi_0
rhoEFlux2
rho_0
p
T
rhoE
Ur2
rhoFlux2
rhoFlux3
thermo:alpha
)
```

`rhoOld`这个变量特别难追踪定位。是在`createDualTimeSteppingFields.H`中定义和赋值，并在`solveUnsteadyFluid.H`中使用的。

最终`forwardStep`算例可以运行，不过运行得很慢。



# 参考文献

1. [DensityBasedTurbo的ppt](http://www.appliedccm.com/wp-content/uploads/2015/07/oliver_borm_slides.pdf)
2. [Implementation of density-based solver for all speeds in the framework of OpenFOAM](http://www.sciencedirect.com/science/article/pii/S0010465514002136)

# 附录

## MRFZone详细对比

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
