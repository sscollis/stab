        module stencils

!=============================================================================!
!.... First derivatives
!=============================================================================!

!.... fourth order central difference ( 1 2 x 4 5 )

        real, parameter :: ga1 =  8.333333333333333333333E-02
        real, parameter :: ga2 = -6.666666666666666666667E-01
        real, parameter :: ga3 =  6.666666666666666666667E-01
        real, parameter :: ga4 = -8.333333333333333333333E-02

!.... fourth order one-sided ( x 2 3 4 5 ) 

        real, parameter :: gc1 = -2.083333333333333333333E+00
        real, parameter :: gc2 =  4.000000000000000000000E+00
        real, parameter :: gc3 = -3.000000000000000000000E+00
        real, parameter :: gc4 =  1.333333333333333333333E+00
        real, parameter :: gc5 = -2.500000000000000000000E-01

!.... fourth order biased difference ( 1 x 2 3 4 5 )

        real, parameter :: gb1 = -2.500000000000000000000E-01
        real, parameter :: gb2 = -8.333333333333333333333E-01
        real, parameter :: gb3 =  1.500000000000000000000E+00
        real, parameter :: gb4 = -5.000000000000000000000E-01
        real, parameter :: gb5 =  8.333333333333333333333E-02

!.... sixth order one-sided ( x 2 3 4 5 6 7 )

        real, parameter :: ge1 = -2.450000000000000000000E+00
        real, parameter :: ge2 =  6.000000000000000000000E+00
        real, parameter :: ge3 = -7.500000000000000000000E+00
        real, parameter :: ge4 =  6.666666666666666666667E+00
        real, parameter :: ge5 = -3.750000000000000000000E+00
        real, parameter :: ge6 =  1.200000000000000000000E+00
        real, parameter :: ge7 = -1.666666666666666666667E-01

!.... sixth order one-pt biased ( 1 x 3 4 5 6 7 ) 

        real, parameter :: gf1 = -1.666666666666666666667E-01
        real, parameter :: gf2 = -1.283333333333333333333E+00
        real, parameter :: gf3 =  2.500000000000000000000E+00 
        real, parameter :: gf4 = -1.666666666666666666667E+00
        real, parameter :: gf5 =  8.333333333333333333333E-01
        real, parameter :: gf6 = -2.500000000000000000000E-01
        real, parameter :: gf7 =  3.333333333333333333333E-02

!==============================================================================
!.... The following stencils are Carpenter's stable and third
!.... order accurate boundary treatment for the explicit fourth
!.... order interior scheme.
!==============================================================================

!.... third order one-sided [Carpenter] ( x 2 3 4 5 6) 
  
        real, parameter :: gg1 = -1.87603205562073772290017534999D+00
        real, parameter :: gg2 =  3.18538322557789286703668141674D+00
        real, parameter :: gg3 = -1.81454567943752757247830550034D+00
        real, parameter :: gg4 =  0.59165824105260274421658150053D+00
        real, parameter :: gg5 = -0.101052068000505624644095417028D+00
        real, parameter :: gg6 =  0.0145883364282753087693133500917D+00

!.... third order biased [Carpenter] ( 1 x 3 4 5 ) 
  
        real, parameter :: gh1 = -0.384254232677925402043233692192D+00
        real, parameter :: gh2 = -0.290638947767348681070279077173D+00
        real, parameter :: gh3 =  0.67176478451531541138011989728D+00
        real, parameter :: gh4 =  0.071081659837399872713651693126D+00
        real, parameter :: gh5 = -0.073630718761724245070378308432D+00
        real, parameter :: gh6 =  0.0056774548542830440901194873933D+00

!.... third order biased [Carpenter] ( 1 2 x 4 5 ) 
  
        real, parameter :: gi1 =  0.182885278686826206583549109115D+00
        real, parameter :: gi2 = -1.08001475417455516433881053795D+00
        real, parameter :: gi3 =  0.65787289649662525818641772732D+00
        real, parameter :: gi4 =  0.1776170486891931456381189546D+00
        real, parameter :: gi5 =  0.076779836395827558602005515075D+00
        real, parameter :: gi6 = -0.0151403060939170046712807681567D+00

!.... third order biased [Carpenter] ( 1 2 3 x 5 ) 
  
        real, parameter :: gj1 = -0.0341837103365291857798396404953D+00
        real, parameter :: gj2 =  0.224829025740103121729831613399D+00
        real, parameter :: gj3 = -0.89081233292845396245426338198D+00
        real, parameter :: gj4 =  0.165299947710035014782196870489D+00
        real, parameter :: gj5 =  0.61343955208752529977826815383D+00
        real, parameter :: gj6 = -0.078572482272680288056193615249D+00

!=============================================================================!
!.... Second derivatives
!=============================================================================!

!.... fourth order central difference

        real, parameter :: da1 = -8.333333333333333333333E-02
        real, parameter :: da2 =  1.333333333333333333333E+00
        real, parameter :: da3 = -2.500000000000000000000E+00
        real, parameter :: da4 =  1.333333333333333333333E+00
        real, parameter :: da5 = -8.333333333333333333333E-02
        
!.... third order biased difference

        real, parameter :: db1 =  9.166666666666666666667E-01
        real, parameter :: db2 = -1.666666666666666666667E+00
        real, parameter :: db3 =  5.000000000000000000000E-01
        real, parameter :: db4 =  3.333333333333333333333E-01
        real, parameter :: db5 = -8.333333333333333333333E-02
        
!.... fourth order one-sided (assumes f'=0) [not smooth]

        real, parameter :: dc1 = -5.763888888888888888889E+00
        real, parameter :: dc2 =  8.000000000000000000000E+00
        real, parameter :: dc3 = -3.000000000000000000000E+00
        real, parameter :: dc4 =  8.888888888888888888889E-01
        real, parameter :: dc5 = -1.250000000000000000000E-01
        
!.... third order one-sided

        real, parameter :: dd1 =  2.916666666666666666667E+00
        real, parameter :: dd2 = -8.666666666666666666667E+00
        real, parameter :: dd3 =  9.500000000000000000000E+00
        real, parameter :: dd4 = -4.666666666666666666667E+00
        real, parameter :: dd5 =  9.166666666666666666667E-01

        end module stencils
