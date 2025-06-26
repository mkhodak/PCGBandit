/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           | Copyright (C) 2016 Ehsan Madadi-Kandjani        |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General macros to create cylinder mesh

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(btQuad, ($1b $2b $2t $1t))
define(topQuad, ($1t $4t $3t $2t))
define(bottomQuad, ($1b $2b $3b $4b))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 1;

// Inner square side half
define(s, 0.02713)

// Inner square side curvature
define(sc, 0.03049)

// cylinder radius
define(r, 48.59e-3)

define(B0, 0.21)

define(mu, 0.75856e-3)
define(sigma, 2.76e6)
define(num_BL, 10)
define(dx_main, 1e-3)
define(BL_ratio, 1.2)
define(dh, calc(1/B0*sqrt(mu/sigma))) // L/Ha
define(dx_wall, calc(dh/num_BL))

define(N_exp, calc(int(log(dx_main/dx_wall)/log(BL_ratio)+1)))
define(L_expansion, calc(dx_wall*(1 - BL_ratio**N_exp)/(1 - BL_ratio)))
define(totalRatio, calc(dx_main/dx_wall))
define(invRatio, calc(1/totalRatio))
define(L_center, calc(2*(r-dh-L_expansion)))
define(N_center, calc(int(L_center/dx_main)))
define(L_outer_main, calc((L_center-2*s)/2))
define(N_outer_main, calc(L_outer_main/dx_main))

define(tw,  8.56e-3)
define(Nw, calc(int(tw/dx_wall)))


// second cylinder radius
define(R, calc(r + tw))



// Base z
define(Zb, calc(-20*r))

// Outlet z
define(Zt, calc(20*r))

// Length of cylinder
define(z, calc(-1*Zb+Zt))

// Number of cells at inner square
define(Ns, calc(int(num_BL+N_exp+N_outer_main)))
// Number of cells between inner square and circle
define(Ni, Ns)

// Number of cells in the cylinder height
define(Nz_reduceBy, 0.1)
define(Nz, calc(int(Nz_reduceBy*z/dx_main)))
define(zLengthMain, 0.95)
define(zLengthBL, 0.05)
define(zCellsMain, 0.9)
define(zCellsBL, 0.1)
define(zRatioMain, 1)
define(zRatioBL, 10)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

define(vert, (x$1$2 y$1$2 $3))
define(evert, (ex$1$2 ey$1$2 $3))

// 45 degree points angle
define(a0, -45)
define(a1, -135)
define(a2, 135)
define(a3, 45)

// Half of 45 degree points angle
define(ea0, 0)
define(ea1, -90)
define(ea2, 180)
define(ea3, 90)

define(ca0, calc(cos((pi/180)*a0)))
define(ca1, calc(cos((pi/180)*a1)))
define(ca2, calc(cos((pi/180)*a2)))
define(ca3, calc(cos((pi/180)*a3)))

define(sa0, calc(sin((pi/180)*a0)))
define(sa1, calc(sin((pi/180)*a1)))
define(sa2, calc(sin((pi/180)*a2)))
define(sa3, calc(sin((pi/180)*a3)))

define(cea0, calc(cos((pi/180)*ea0)))
define(cea1, calc(cos((pi/180)*ea1)))
define(cea2, calc(cos((pi/180)*ea2)))
define(cea3, calc(cos((pi/180)*ea3)))

define(sea0, calc(sin((pi/180)*ea0)))
define(sea1, calc(sin((pi/180)*ea1)))
define(sea2, calc(sin((pi/180)*ea2)))
define(sea3, calc(sin((pi/180)*ea3)))

// Inner square x and y position

// x
define(x00, s)
define(x01, calc(-1.0*s))
define(x02, calc(-1.0*s))
define(x03, s)

// y
define(y00, calc(-1.0*s))
define(y01, calc(-1.0*s))
define(y02, s)
define(y03, s)

// Circle x and y positions

// x
define(x10, calc(r*ca0))
define(x11, calc(r*ca1))
define(x12, calc(r*ca2))
define(x13, calc(r*ca3))

// y
define(y10, calc(r*sa0))
define(y11, calc(r*sa1))
define(y12, calc(r*sa2))
define(y13, calc(r*sa3))


// Second Circle x and y positions

// x
define(x20, calc(R*ca0))
define(x21, calc(R*ca1))
define(x22, calc(R*ca2))
define(x23, calc(R*ca3))

// y
define(y20, calc(R*sa0))
define(y21, calc(R*sa1))
define(y22, calc(R*sa2))
define(y23, calc(R*sa3))

// Inner square x and y position middle curvatures

// x
define(ex00, sc)
define(ex01, 0)
define(ex02, calc(-1.0*sc))
define(ex03, 0)

// y
define(ey00, 0)
define(ey01, calc(-1.0*sc))
define(ey02, 0)
define(ey03, sc)

// Circle x and y positions middle curvatures

// x
define(ex10, calc(r*cea0))
define(ex11, calc(r*cea1))
define(ex12, calc(r*cea2))
define(ex13, calc(r*cea3))

// y
define(ey10, calc(r*sea0))
define(ey11, calc(r*sea1))
define(ey12, calc(r*sea2))
define(ey13, calc(r*sea3))


// Second Circle x and y positions middle curvatures

// x
define(ex20, calc(R*cea0))
define(ex21, calc(R*cea1))
define(ex22, calc(R*cea2))
define(ex23, calc(R*cea3))

// y
define(ey20, calc(R*sea0))
define(ey21, calc(R*sea1))
define(ey22, calc(R*sea2))
define(ey23, calc(R*sea3))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    vert(0, 0, Zb) vlabel(s0b)
    vert(0, 1, Zb) vlabel(s1b)
    vert(0, 2, Zb) vlabel(s2b)
    vert(0, 3, Zb) vlabel(s3b)
    
    vert(1, 0, Zb) vlabel(r0b)
    vert(1, 1, Zb) vlabel(r1b)
    vert(1, 2, Zb) vlabel(r2b)
    vert(1, 3, Zb) vlabel(r3b)
    
    vert(0, 0, Zt) vlabel(s0t)
    vert(0, 1, Zt) vlabel(s1t)
    vert(0, 2, Zt) vlabel(s2t)
    vert(0, 3, Zt) vlabel(s3t)
    
    vert(1, 0, Zt) vlabel(r0t)
    vert(1, 1, Zt) vlabel(r1t)
    vert(1, 2, Zt) vlabel(r2t)
    vert(1, 3, Zt) vlabel(r3t)

    vert(2, 0, Zb) vlabel(w0b)
    vert(2, 1, Zb) vlabel(w1b)
    vert(2, 2, Zb) vlabel(w2b)
    vert(2, 3, Zb) vlabel(w3b)

    vert(2, 0, Zt) vlabel(w0t)
    vert(2, 1, Zt) vlabel(w1t)
    vert(2, 2, Zt) vlabel(w2t)
    vert(2, 3, Zt) vlabel(w3t)


);

blocks
(
    //block0
    hex2D(s1, s0, s3, s2)
    liquid
    (Ns Ns Nz)
    simpleGrading (1 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )
    
    //block1
    hex2D(s0, r0, r3, s3)
    liquid
    (Ni Ns Nz)
    simpleGrading (((L_outer_main N_outer_main 1)(L_expansion N_exp invRatio)(dh num_BL 1)) 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )
    
    //block2
    hex2D(s3, r3, r2, s2)
    liquid
    (Ni Ns Nz)
    simpleGrading (((L_outer_main N_outer_main 1)(L_expansion N_exp invRatio)(dh num_BL 1)) 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )
    
    //block3
    hex2D(s2, r2, r1, s1)
    liquid
    (Ni Ns Nz)
    simpleGrading (((L_outer_main N_outer_main 1)(L_expansion N_exp invRatio)(dh num_BL 1)) 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )
    
    //block4
    hex2D(s1, r1, r0, s0)
    liquid
    (Ni Ns Nz)
    simpleGrading (((L_outer_main N_outer_main 1)(L_expansion N_exp invRatio)(dh num_BL 1)) 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )

    //block5
    hex2D(r0, w0, w3, r3)
    solidWalls
    (Nz Ns Nz)
    simpleGrading (1 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )

    //block6
    hex2D(r1, w1, w0, r0)
    solidWalls
    (Nz Ns Nz)
    simpleGrading (1 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )

    //block7
    hex2D(r2, w2, w1, r1)
    solidWalls
    (Nz Ns Nz)
    simpleGrading (1 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )

    //block8
    hex2D(r3, w3, w2, r2)
    solidWalls
    (Nz Ns Nz)
    simpleGrading (1 1 ((zLengthBL zCellsBL zRatioBL)(zLengthMain zCellsMain zRatioMain)) )


);

edges
(
    //Circle edges
    arc r3b r0b evert(1, 0, Zb)
    arc r0b r1b evert(1, 1, Zb)
    arc r1b r2b evert(1, 2, Zb)
    arc r2b r3b evert(1, 3, Zb)
    
    //Circle edges
    arc r3t r0t evert(1, 0, Zt)
    arc r0t r1t evert(1, 1, Zt)
    arc r1t r2t evert(1, 2, Zt)
    arc r2t r3t evert(1, 3, Zt)
    
    arc s3b s0b evert(0, 0, Zb)
    arc s0b s1b evert(0, 1, Zb)
    arc s1b s2b evert(0, 2, Zb)
    arc s2b s3b evert(0, 3, Zb)
    
    arc s3t s0t evert(0, 0, Zt)
    arc s0t s1t evert(0, 1, Zt)
    arc s1t s2t evert(0, 2, Zt)
    arc s2t s3t evert(0, 3, Zt)
    

    //Second Circle edges
    arc w3b w0b evert(2, 0, Zb)
    arc w0b w1b evert(2, 1, Zb)
    arc w1b w2b evert(2, 2, Zb)
    arc w2b w3b evert(2, 3, Zb)
    
    arc w3t w0t evert(2, 0, Zt)
    arc w0t w1t evert(2, 1, Zt)
    arc w1t w2t evert(2, 2, Zt)
    arc w2t w3t evert(2, 3, Zt)



);


boundary
(
    exteriorWalls
    {
    type patch;
    faces
    (
        btQuad(w0, w3)
        btQuad(w1, w0)
        btQuad(w2, w1)
        btQuad(w3, w2)
        bottomQuad(r3, w3, r2, w2)
        bottomQuad(r2, w2, r1, w1)
        bottomQuad(r1, w1, r0, w0)
        bottomQuad(r0, w0, r3, w3)
        topQuad(r3, w3, r2, w2)
        topQuad(r2, w2, r1, w1)
        topQuad(r1, w1, r0, w0)
        topQuad(r0, w0, r3, w3)
    );
    }

    inlet
    {
    type patch;
    faces 
    (
        bottomQuad(s3, s0, s1, s2)
        bottomQuad(s3, r3, r0, s0)
        bottomQuad(s2, r2, r3, s3)
        bottomQuad(s1, r1, r2, s2)
        bottomQuad(s0, r0, r1, s1)
    );
    }
    
    sink
    {
    type patch;
    faces 
    (
        topQuad(s3, s0, s1, s2)
        topQuad(s3, r3, r0, s0)
        topQuad(s2, r2, r3, s3)
        topQuad(s1, r1, r2, s2)
        topQuad(s0, r0, r1, s1)
    );
    }
);

mergePatchPairs
(
);
