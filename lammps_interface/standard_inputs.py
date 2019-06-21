
input_files = {

'reaxff': """    units        real
    atom_style   full

    read_data "lmp.data"


    pair_style reax/c control.reaxff
    pair_coeff * *  {}  O H Ti

    fix             2 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

    minimize 1.0e-5 1.0e-7 1000 10000
    timestep {}
    compute energy all pe/atom
    variable energy equal c_energy
    dump molfile all custom 500 atoms.atm type x y z fx fy fz c_energy
    #fix 1 all 1 0.0 10.0 1.0e-6 qeq/reax
    #fix   fxnvt all nvt temp 300.0 300.0 500.0 tchain 1
    fix   fxnvp all npt temp {} 300.0 500.0 tchain 1 press 1.0 1.0 1
    #compute myRDF all rdf 100
    #fix 1 all ave/time 100 1 1000000 c_myRDF[*] file tmp.rdf mode vector
    run   {}""",

'reaxff_diffusion': """    units        real
    atom_style   full

    read_data "mt.data"


    pair_style reax/c control.reaxff
    pair_coeff * *  2_ffield.reax.water_2017  O H

    fix             2 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c
    compute energy all pe/atom
    variable energy equal c_energy

    dump molfile all custom 500 atoms.atm type x y z fx fy fz c_energy

    minimize 1.0e-5 1.0e-7 1000 10000
    timestep {}

    fix             4 all npt temp 300.0 300.0 100.0 iso 0.1 1 10

    run 5000

    unfix 4

    #reset_timestep 0

    compute energy all pe/atom
    variable energy equal c_energy
    #dump molfile all custom 500 atoms.atm type x y z fx fy fz c_energy
    #fix 1 all 1 0.0 10.0 1.0e-6 qeq/reax
    fix   fxnvt all nvt temp {} 300.0 100.0 tchain 1
    compute         msd all msd com yes
    variable        twopoint equal c_msd[4]/6/(step*dt+1.0e-6)
    fix             9 all vector 10 c_msd[4]
    variable        fitslope equal slope(f_9)/6/(10*dt)

    thermo_style    custom step temp c_msd[4] v_twopoint v_fitslope

    # only need to run for 10K steps to make a good 100-frame movie

    #dump           1 all custom 1 tmp.dump id type vx vy vz

    #dump           2 all image 100 image.*.jpg type type zoom 1.6 adiam 1.2

    thermo          1000
    run   {}""",


}
