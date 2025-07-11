get_view(ls::LScene) = ls.scene.camera.view[]
set_view(ls::LScene, new::AbstractMatrix) = (ls.scene.camera.view[] = new)

function set_orthographic(ax::LScene)
    cam = cameracontrols(ax.scene)
    cam.fov[] = 1
    return nothing
end

hide_axis(ls::LScene) = (ls.show_axis[] = false)

function take_screenshot(
    fname::String,
    system::BMO.Nullable{BMO.AbstractSystem},
    beam::BMO.Nullable{BMO.AbstractBeam};
    size::Tuple{Int, Int} = (600, 300),
    view::BMO.Nullable{AbstractMatrix} = nothing,
    flen::Real = 10e-2,
    px_per_unit::Int = 8,
    color::Union{Symbol, RGBf, RGBAf} = :red,
    optional_rendering_fct::BMO.Nullable{Function} = nothing
    )
    fig = Figure(; size)
    display(fig)
    ax = LScene(fig[1,1])
    if !isnothing(system)
        render!(ax, system)
    end
    if !isnothing(beam)
        render!(ax, beam; flen, color)
    end
    if !isnothing(view)
        set_view(ax, view)
    end
    if !isnothing(optional_rendering_fct)
        optional_rendering_fct(fig, ax)
    end
    hide_axis(ax)
    save(fname, fig; px_per_unit, update = false)
end

arrow!(
    ax::LScene,
    pos::AbstractVector,
    dir::AbstractVector;
    # Custom kwargs
    scale=1,
    # Makie kwargs
    color=:blue,
    tiplength=0.2,
    tipradius=0.1,
    ) = arrows3d!(
        ax,
        [Point3(pos)],
        [Point3(dir*5e-3*scale)];
        tiplength,
        tipradius,
        color
)

lens_color() = RGBf(0.678, 0.847, 0.902)
lens_color(alpha) = RGBAf(0.678, 0.847, 0.902, alpha)