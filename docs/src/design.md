# API design

!!! warning
    This page is very much WIP...

## Conventions

### Right-handedness

All **coordinate systems** are or must be defined **right-handed**! All normal vectors are or must be defined right-handed! 
All **rotations** are or must be performed in a **counter-clockwise** manner for a positive rotation angle ``\theta > 0`` and vice-versa!
For a definition of rotation matrix order, refer to this [article](https://dominicplein.medium.com/extrinsic-intrinsic-rotation-do-i-multiply-from-right-or-left-357c38c1abfd).

!!! warning
    Failure to comply with this convention can lead to spurious effects and silent bugs when using the kinematic API of this package!

## Intersect - Interact - Repeat - Loop (IIRP)

### Tracing

### Retracing

### CPU and GPU

!!! info
    GPU processing (tracing) of optical systems is  not supported at the moment.

## Geometry representation

### Meshes

### Signed Distance Functions (SDFs)

For an introduction into SDFs the [website of Inigo Quilez](https://iquilezles.org/articles/distfunctions/) is referred to. The following shapes have been implemented:

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractSDF);
```