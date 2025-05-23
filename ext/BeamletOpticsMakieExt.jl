module BeamletOpticsMakieExt

using BeamletOptics: faces, vertices, AbstractMesh, AbstractRay, Beam, PolarizedRay, intensity, rays, NonInteractableObject
import BeamletOptics: render_object!, render_ray!, _render_beam!,
    _render_ray!, render_surface!, _render_object_normal!, render_sdf_mesh!, render_dummy_mesh!
using Makie: Axis3, LScene, mesh!, surface!, lines!, RGBAf, scatter!
using GeometryBasics: Point2, Point3
using AbstractTrees: PreOrderDFS

const _RenderEnv = Union{Axis3, LScene}

function render_object!(axis::_RenderEnv, mesh::AbstractMesh; transparency=true)
    mesh!(axis, vertices(mesh), faces(mesh); transparency)
    return nothing
end

function _render_ray!(axis::_RenderEnv,
        ray::AbstractRay,
        ray_end::AbstractVector;
        color = :blue,
        show_pos = false)
    lines!(axis,
        [ray.pos[1], ray_end[1]],
        [ray.pos[2], ray_end[2]],
        [ray.pos[3], ray_end[3]],
        color=color,
        linewidth=1.0,
        transparency=true)
    if show_pos
        scatter!(axis, ray.pos; color)
    end
    return nothing
end

function _render_beam!(axis, beam::Beam{T, R}; color=:red, flen=1.0) where {T<:Real, R<:PolarizedRay{T}}
    I0 = intensity(first(rays(beam)))
    for child in PreOrderDFS(beam)
        for ray in rays(child)
            I = intensity(ray)
            render_ray!(axis, ray, color = RGBAf(1,0,0, I/I0), flen = flen)
        end
    end
    return nothing
end

function render_surface!(axis::_RenderEnv, X, Y, Z; kwargs...)
    surface!(axis, X, Y, Z; kwargs...)
end

function _render_object_normal!(axis::_RenderEnv,
        pos::AbstractVector,
        vec::AbstractVector;
        color = :blue)
    lines!(axis,
        [pos[1], vec[1]],
        [pos[2], vec[2]],
        [pos[3], vec[3]],
        color=color)
end

render_sdf_mesh!(axis::_RenderEnv, vertices, faces; transparency = true) = mesh!(axis, vertices, faces, transparency=transparency)

function render_dummy_mesh!(axis::_RenderEnv, d::NonInteractableObject; transparency = false)
    mesh = d.shape
    mesh!(axis, vertices(mesh), faces(mesh); transparency, color = :grey)
    return nothing
end

end
