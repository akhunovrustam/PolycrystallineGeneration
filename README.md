# README #

This is repo for master thesis of Akhunov Rustam. Supervisor: Szymon Winczewski.

The topic of the thesis is building of polycrystal lattice from the distribution of sizes and orientation of grains.

In order to build the project you need to install Voro++ library. Then you can from the root of project run:
	
	g++ face_to_place.cc -L/usr/local/lib -lvoro++ -I/usr/local/include/voro++

Note that command works only on Posix systems.
After running the binary you will obtain 4 files. From 2 of them you can see grains of Voronoi tessellation by typing: 

	splot "random_points_p.gnu" u 2:3:4, "random_points_v.gnu" with lines
	
	
Other 2 files present lattice itself and points of the vertices of the tessellation in XYZ format. 
You can display them in VMD for example