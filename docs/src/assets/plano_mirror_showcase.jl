using CairoMakie, BeamletOptics

# relative to the .md file location from which this code is included
asset_dir = joinpath(@__DIR__, "..", "assets")

function spawn_mirror_mount()
    m1 = RoundPlanoMirror(1BeamletOptics.inch, 6e-3)
    holder = MeshDummy(joinpath(asset_dir, "Mirror_Post.stl"))
    translate_to3d!(holder, [0,0,-(5.68e-2)])
    BeamletOptics.set_new_origin3d!(holder)
    mirror_mount = ObjectGroup([holder, m1])
    return mirror_mount
end

m1 = spawn_mirror_mount()
m2 = spawn_mirror_mount()
m3 = spawn_mirror_mount()
m4 = spawn_mirror_mount()

zrotate3d!(m1, deg2rad(45))
zrotate3d!(m2, deg2rad(-135))
zrotate3d!(m3, deg2rad(-45))
zrotate3d!(m4, deg2rad(135))

dy = 0.3

translate_to3d!(m1, [0,    0,      0])
translate_to3d!(m2, [0.1,  0,      0])
translate_to3d!(m3, [0.1,  dy,     0])
translate_to3d!(m4, [0,    dy,     0])

system = StaticSystem([m1, m2, m3, m4])

fig = Figure(size=(600, 300))
limits = (-0.05, 0.15, -0.1, 0.4, -0.1, 0.05)
aspect = (2,5,1.5)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0.0, elevation=pi/2)

hidedecorations!(ax)
hidespines!(ax)

beam = GaussianBeamlet([0, -0.2, 0], [0, 1, 0], 532e-9, 1e-3)

solve_system!(system, beam, r_max=100)

render_object!(ax, m1)
render_object!(ax, m4)
render_object!(ax, m2)
render_object!(ax, m3)
render_beam!(ax, beam, flen=0.2, color=RGBf(0,1,0))
