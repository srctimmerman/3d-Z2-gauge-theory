#This is the data structure that relates up-cube to spin indices,
#   i.e. given an index labeling an up-cube, this holds all of the other
#   spin indices in that up-cube, organized into xy, xz, and yz plaquettes.
#   Note that only indices on the x-lattice label up-cubes.
function make_up_cube(L, U=L^2, V=L^3)
    UpCube = zeros(Int,V,3,4)
    for i = 1:V

        #neighboring spin indices in the up cube
        x = i            #spins living on the x lattice
        xy = x + L
        xz = x + U

        y = i + V        #spins living on the y lattice
        yx = y + 1
        yz = y + U

        z = i + 2*V      #spins living on the z lattice
        zx = z + 1
        zy = z + L

        #enforce periodic boundary conditions
        if mod(i,L) == 0         #fix x-direction
            yx -= L
            zx -=L
        end
        if mod(i-1, U) + L >= U  #fix y-direction
            xy -= U
            zy -= U
        end
        if i + U > V             #fix z-direction
            xz -= V
            yz -= V
        end

        #fill up UpCube
        UpCube[i,1,1] = x
        UpCube[i,1,2] = y  #xy plaquette
        UpCube[i,1,3] = xy
        UpCube[i,1,4] = yx

        UpCube[i,2,1] = x
        UpCube[i,2,2] = z  #xz plaquette
        UpCube[i,2,3] = xz
        UpCube[i,2,4] = zx

        UpCube[i,3,1] = y
        UpCube[i,3,2] = zy  #yz plaquette
        UpCube[i,3,3] = yz
        UpCube[i,3,4] = z
    end
    return UpCube
end

# #for debugging
# L = 3
# UpCube = make_up_cube(L)
# for i = 1:L^3
#     print(i," ", UpCube[i,3,1], " ", UpCube[i,3,2], " ", UpCube[i,3,3], " ", UpCube[i,3,4], "\n")
# end


#This is the inverse data structure that relates a spin index to its 4 plaquettes
#   i.e. given a spin index, this finds all of the plaquettes that involve that spin
function make_associated_cube(L)

    N = 3*L^3
    AssociatedCube = zeros(Int,N,4,4)

    for i = 1:N
        #counter for how many plaquettes we've found for this spin site
        p = 1

        #iterate over UpCube to find the plaquettes involving i
        for j = 1:L^3
            for k = 1:3
                for l = 1:4
                    if UpCube[j,k,l] == i    #we've found a plaquette!
                        AssociatedCube[i,p,:] = UpCube[j,k,:]
                        p+=1
                    end
                end
            end
        end

    end
    return AssociatedCube
end

#for debugging
# L = 3
# AssociatedCube = make_associated_cube(L);
# for i = 1:4
#     print(i," ", AssociatedCube[1,i,1], " ", AssociatedCube[1,i,2], " ", AssociatedCube[1,i,3], " ", AssociatedCube[1,i,4], "\n")
# end

# this is a data structure that, given a LATTICE index (a vertex), finds all of
#   the spins living on the bonds connected to that vertex (i.e. the star)
function make_cube_stars(L, U=L^2, V=L^3)
        CubeStars = zeros(Int,V,6)
        for i = 1:V

            xp = i
            xm = xp - 1
            yp = i + V
            ym = yp - L
            zp = i + 2*V
            zm = zp - L^2

            #enforce periodic boundary conditions
            if mod(i,L) == 1         #fix x-direction: i was in first column
                xm += L
            end
            if mod(i-1, U) < L  ``     #fix y-direction: i was in first row
                ym += L^2
            end
            if i <= L^2             #fix z-direction: i was in first layer
                zm += L^3
            end

            #fill the cube star
            CubeStars[i,1] = xp
            CubeStars[i,2] = xm
            CubeStars[i,3] = yp
            CubeStars[i,4] = ym
            CubeStars[i,5] = zp
            CubeStars[i,6] = zm
        end
        return CubeStars
end

# # for debugging
# L = 3
# CubeStars = make_cube_stars(L);
# for i = 1:L^3
#     print(CubeStars[i,1]," ", CubeStars[i,2], " ", CubeStars[i,3], " ", CubeStars[i,4], " ", CubeStars[i,5], " ", CubeStars[i,6], "\n")
# end
