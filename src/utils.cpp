#include "utils.hpp"



float3x3 rotation_matrix(float3 rotation_angle){
    float3x3 rotation_matrix;
    float3x3 rotation_x = float3x3(1, 0, 0, 
    0, cos(rotation_angle.x), -sin(rotation_angle.x), 
    0, sin(rotation_angle.x), cos(rotation_angle.x));
    float3x3 rotation_y = float3x3(cos(rotation_angle.y), 0, sin(rotation_angle.y),
    0, 1, 0,
    -sin(rotation_angle.y), 0, cos(rotation_angle.y));
    float3x3 rotation_z = float3x3(cos(rotation_angle.z), -sin(rotation_angle.z), 0,
    sin(rotation_angle.z), cos(rotation_angle.z), 0,
    0, 0, 1);
    rotation_matrix = rotation_x * rotation_y * rotation_z;
    return rotation_matrix;
}

string id2object(int id, float3 rotation_angle, LBM & lbm, float size) {
    string shape_name;
    float3 obj_center = lbm.center();
    float3x3 rotation = rotation_matrix(rotation_angle);
    
    
    // set the shape name based on the id
    switch(id) {
        case 0:
            shape_name = "dog";
            lbm.voxelize_stl(get_exe_path()+"../stl/Dog_Singlecolor.stl", obj_center, rotation, size);
            break;
        case 1:
            shape_name = "shark";
            lbm.voxelize_stl(get_exe_path()+"../stl/jeff_the_land_shark.stl", obj_center, rotation, size);
            break;
        case 2:
            shape_name = "landing gear";
            lbm.voxelize_stl(get_exe_path()+"../stl/crm-hl_reference_ldg.stl", obj_center, rotation, size);
            break;
        case 3:
            shape_name = "cat ring";
            lbm.voxelize_stl(get_exe_path()+"../stl/cat ring.stl", obj_center, rotation, size);
            break;
        case 4:
            shape_name = "flexible cat";
            lbm.voxelize_stl(get_exe_path()+"../stl/CatFlexi.stl", obj_center, rotation, size);
            break;
        case 5:
            shape_name = "gong-gi";
            lbm.voxelize_stl(get_exe_path()+"../stl/gong-gi.stl", obj_center, rotation, size);
            break;
        case 6:
            shape_name = "octopus";
            lbm.voxelize_stl(get_exe_path()+"../stl/Octopus_spiral_v6.stl", obj_center, rotation, size);
            break;
        case 7:
            shape_name = "planter";
            lbm.voxelize_stl(get_exe_path()+"../stl/radiant_planter.stl", obj_center, rotation, size);
            break;
        case 8:
            shape_name = "plate";
            lbm.voxelize_stl(get_exe_path()+"../stl/radiant_plate.stl", obj_center, rotation, size);
            break;
        case 9:
            shape_name = "shark v2.0";
            lbm.voxelize_stl(get_exe_path()+"../stl/Shark_v2.0_A.stl", obj_center, rotation, size);
            break;
        case 10:
            shape_name = "car wheel SKODA";
            lbm.voxelize_stl(get_exe_path()+"../stl/SKODA_ENYAQ_FL_COUPE_S24_WHEEL_C5W_R21_20241220.stl", obj_center, rotation, size);
            break;
        case 11:
            shape_name = "Zorro character";
            lbm.voxelize_stl(get_exe_path()+"../stl/zorro_chibi2fix.stl", obj_center, rotation, size);
            break;
        case 12:
            shape_name = "Cow";
            lbm.voxelize_stl(get_exe_path()+"../stl/Cow_t.stl", obj_center, rotation, size);
            break;
        case 13:
            shape_name = "cone_0";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cone_0.stl", obj_center, rotation, size);
            break;
        case 14:
            shape_name = "cone_1";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cone_1.stl", obj_center, rotation, size);
            break;
        case 15:
            shape_name = "cone_2";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cone_2.stl", obj_center, rotation, size);
            break;
        case 16:    
            shape_name = "cone_3";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cone_3.stl", obj_center, rotation, size);
            break;
        case 17:
            shape_name = "cube_0";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cube_0.stl", obj_center, rotation, size);
            break;
        case 18:
            shape_name = "cube_1";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cube_1.stl", obj_center, rotation, size);
            break;
        case 19:
            shape_name = "cube_2";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cube_2.stl", obj_center, rotation, size);
            break;
        case 20:
            shape_name = "cube_3";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cube_3.stl", obj_center, rotation, size);
            break;
        case 21:
            shape_name = "cylinder_0";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cylinder_0.stl", obj_center, rotation, size);
            break;
        case 22:
            shape_name = "cylinder_1";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cylinder_1.stl", obj_center, rotation, size);
            break;
        case 23:    
            shape_name = "cylinder_2";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cylinder_2.stl", obj_center, rotation, size);
            break;
        case 24:
            shape_name = "cylinder_3";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/cylinder_3.stl", obj_center, rotation, size);
            break;
        case 25:    
            shape_name = "ellipsoid_0";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/sphere_0.stl", obj_center, rotation, size);
            break;
        case 26:
            shape_name = "ellipsoid_1";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/sphere_1.stl", obj_center, rotation, size);
            break;
        case 27:
            shape_name = "ellipsoid_2";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/sphere_2.stl", obj_center, rotation, size);
            break;
        case 28:
            shape_name = "rectangular_0";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/rectangular_0.stl", obj_center, rotation, size);
            break;
        case 29:
            shape_name = "rectangular_1";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/rectangular_1.stl", obj_center, rotation, size);
            break;
        case 30:
            shape_name = "rectangular_2";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/rectangular_2.stl", obj_center, rotation, size);
            break;
        case 31:
            shape_name = "sphere_0";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/sphere_0.stl", obj_center, rotation, size);
            break;
        case 32:
            shape_name = "sphere_1";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/sphere_1.stl", obj_center, rotation, size);
            break;
        case 33:
            shape_name = "sphere_2";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/sphere_2.stl", obj_center, rotation, size);
            break;
        case 34:
            shape_name = "sphere_3";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/sphere_3.stl", obj_center, rotation, size);
            break;
        case 35:
            shape_name = "torus_0";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/torus_0.stl", obj_center, rotation, size);
            break;
        case 36:
            shape_name = "torus_1";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/torus_1.stl", obj_center, rotation, size);
            break;
        case 37:
            shape_name = "torus_2";
            lbm.voxelize_stl(get_exe_path()+"../generated_stl/torus_2.stl", obj_center, rotation, size);
            break;
        default:
            shape_name = "unknown object";
            break;
    }

	// construct the prompt
    string twist_info = std::to_string(int(rotation_angle.x*100)) + ", " + std::to_string(int(rotation_angle.y*100)) + ", " + std::to_string(int(rotation_angle.z*100));
	string prompt = "(" + twist_info + "; " + shape_name + ")" + " Visualize the Q-criterion of the " + shape_name + " to represent the turbulent flow field of the object in the given image.";

	return prompt;
}

void create_dataset(uint frequency,bool is_train, LBM &lbm, string prompt){
    string original_path = "./";
    // string original_path = "C:\\Users\\86177\\Desktop\\workspace\\EXP\\data_generation\\";

    const uint star_T = 10000u; // number of LBM time steps to simulate
	const uint frames = 49; // number of LBM time steps to simulate
	const uint lbm_T = star_T + frames*frequency; // number of LBM time steps to simulate
	
    size_t start_pos = prompt.find('('); // 找到第一个 '(' 的位置
    size_t end_pos = prompt.find(')');   // 找到第一个 ')' 的位置
    uint Nx = lbm.get_Nx();
    uint Ny = lbm.get_Ny();
    uint Nz = lbm.get_Nz();
    string object_id_str = prompt.substr(start_pos + 1, end_pos - start_pos - 1); // 提取括号内的内容
    string folder;
    if (is_train){
        folder = "train";
    }else{
        folder = "test";
    }
    // save the first frame
    lbm.run(0, lbm_T); // initialize simulation
    lbm.graphics.set_camera_free(
    float3(((float)Nz/Nx) * Nx * 1.9f, 0.0f * Ny, 0.0f * Nz),  // 相机位置（x, y, z）
    0.0f,   // 相机朝向的偏航角度 (yaw)
    0.0f,   // 相机朝向的俯仰角度 (pitch)
    50.0f   // 相机视距
    );
    lbm.graphics.write_frame(original_path + "SynData/" + folder + "/images/" + object_id_str + "-"); // export image from camera position
    lbm.run(star_T, lbm_T); // initialize simulation

	while(lbm.get_t()<lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 25.0f)) { // render enough frames for 25 seconds of 60fps video
			// 设置相机位置在管道侧边，并正对管道
			lbm.graphics.set_camera_free(
				float3(((float)Nz/Nx) * Nx * 1.9f, 0.0f * Ny, 0.0f * Nz), // 相机位置（x, y, z）
				0.0f,   // 相机朝向的偏航角度 (yaw)
				0.0f,   // 相机朝向的俯仰角度 (pitch)
				50.0f   // 相机视距
			);
			// 生成完整的路径
			std::string camera_path = original_path + "SynData/" + folder + "/images/" + object_id_str + "/";
			// MAKE_DIR(camera_path);
			lbm.graphics.write_frame(camera_path); // export image from camera position 
		}
		lbm.run(frequency, lbm_T); // run 1 LBM time step
	}
	// open the SynData/" + folder + "/videos.txt/prompts.txt write the prompt
	std::ofstream file(original_path + "SynData/" + folder + "/prompts.txt", std::ios::app);
	file <<  prompt + "\n";
	file.close();
	// open the SynData/" + folder + "/videos.txt/videos.txt write the video path
	std::ofstream file1(original_path + "SynData/" + folder + "/videos.txt", std::ios::app);
	file1 <<"videos/" + object_id_str + ".mp4\n";
	file1.close();

    // open the SynData/" + folder + "/videos.txt/images.txt write the video path
    std::ofstream file2(original_path + "SynData/" + folder + "/images.txt", std::ios::app);
    file2 << "images/" + object_id_str + "-" + "image-000000000.png\n";
    file2.close();
}


void run_simulation(float si_u, int id, float3 rotation, float size, bool is_train, bool render = true) {
    const uint3 lbm_N = resolution(float3(0.01f, 2.0f, 1.0f), 1000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	const float si_length = 2.4f;
	// const float si_T = 10.0f;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	const float lbm_length = 0.65f*(float)lbm_N.y; // length of the simulation box in lbm
	const float lbm_u = 0.075f; // velocity of the fluid in lbm
	units.set_m_kg_s(lbm_length, lbm_u, 1.0f, si_length, si_u, si_rho);        
	const float lbm_nu = units.nu(si_nu);
	// const ulong lbm_T = units.t(si_T);
	print_info("Re = "+std::to_string(to_uint(units.si_Re(si_length, si_u, si_nu))));
    // initialize the LBM object
	LBM lbm(lbm_N, lbm_nu);
    const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();
    size = size * Nz;
    string prompt = id2object(id, rotation, lbm, size);

    // find the boundary location
    const uint Nx_1 = Nx - 1;
    const uint Ny_1 = Ny - 1;
    const uint Nz_1 = Nz - 1;
    parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		//if(z==0u) lbm.flags[n] = TYPE_S; // solid floor
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u; // initialize y-velocity everywhere except in solid cells
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u||z==0u) lbm.flags[n] = TYPE_E; // all other simulation box boundaries are inflow/outflow
	}); // ####################################################################### run simulation, export images and data ##########################################################################
    // set the Visualize flag
    lbm.graphics.visualization_modes = VIS_Q_CRITERION|VIS_FLAG_SURFACE;
    if (render){

        float dt = units.si_t(1ul); // 1 LBM步对应的真实时间，单位是秒
        print_info("dt = " + std::to_string(dt));
        float slow_down = 100.0f; 
        int frequency = (int)((1/24.0)/dt/slow_down); // 24fps
        // int frequency = 1;
        print_info("frequency = " + std::to_string(frequency));
        create_dataset(frequency, is_train, lbm, prompt);
    }else{
        lbm.run();
    }

}

void create_dataset_2d(uint frequency,bool is_train, LBM &lbm, string prompt){
    string original_path = "./";
    // string original_path = "C:\\Users\\86177\\Desktop\\workspace\\EXP\\data_generation\\";

    const uint star_T = 100000u; // number of LBM time steps to simulate
	const uint frames = 49; // number of LBM time steps to simulate
	const uint lbm_T = star_T + frames*frequency; // number of LBM time steps to simulate
	
    size_t start_pos = prompt.find('('); // 找到第一个 '(' 的位置
    size_t end_pos = prompt.find(')');   // 找到第一个 ')' 的位置
    uint Nx = lbm.get_Nx();
    uint Ny = lbm.get_Ny();
    uint Nz = lbm.get_Nz();
    string object_id_str = prompt.substr(start_pos + 1, end_pos - start_pos - 1); // 提取括号内的内容
    string folder;
    if (is_train){
        folder = "train";
    }else{
        folder = "test";
    }
    // save the first frame
    lbm.run(0, lbm_T); // initialize simulation
    lbm.graphics.set_camera_free(
    float3(0.0f * Nx, 0.0f * Ny, ((float)Nx/Nz) * Nz * 1.8f),  // 相机位置（x, y, z）
    0.0f,   // 相机朝向的偏航角度 (yaw)
    0.0f,   // 相机朝向的俯仰角度 (pitch)
    50.0f   // 相机视距
    );
    lbm.graphics.write_frame(original_path + "SynData_2D/" + folder + "/images/" + object_id_str + "-"); // export image from camera position
    lbm.run(star_T, lbm_T); // initialize simulation

	while(lbm.get_t()<lbm_T) { // main simulation loop
		if(lbm.graphics.next_frame(lbm_T, 25.0f)) { // render enough frames for 25 seconds of 60fps video
			// 设置相机位置在管道侧边，并正对管道
			lbm.graphics.set_camera_free(
				float3(0.0f * Nx, 0.0f * Ny, ((float)Nx/Nz) * Nz * 1.8f),  // 相机位置（x, y, z）
				0.0f,   // 相机朝向的偏航角度 (yaw)
				0.0f,   // 相机朝向的俯仰角度 (pitch)
				50.0f   // 相机视距
			);
			// 生成完整的路径
			std::string camera_path = original_path + "SynData_2D/" + folder + "/images/" + object_id_str + "/";
			// MAKE_DIR(camera_path);
			lbm.graphics.write_frame(camera_path); // export image from camera position 
		}
		lbm.run(frequency, lbm_T); // run 1 LBM time step
	}
	// open the SynData/" + folder + "/videos.txt/prompts.txt write the prompt
	std::ofstream file(original_path + "SynData_2D/" + folder + "/prompts.txt", std::ios::app);
	file <<  prompt + "\n";
	file.close();
	// open the SynData/" + folder + "/videos.txt/videos.txt write the video path
	std::ofstream file1(original_path + "SynData_2D/" + folder + "/videos.txt", std::ios::app);
	file1 <<"videos/" + object_id_str + ".mp4\n";
	file1.close();

    // open the SynData/" + folder + "/videos.txt/images.txt write the video path
    std::ofstream file2(original_path + "SynData_2D/" + folder + "/images.txt", std::ios::app);
    file2 << "images/" + object_id_str + "-" + "image-000000000.png\n";
    file2.close();
}

bool id2object_2d(int id, int x, int y, int z, LBM &lbm, float size, string &prompt) {
    float3 obj_center = float3(lbm.get_Nx()/2, lbm.get_Ny()/2, lbm.get_Nz()/2);
    switch(id) {
        case 0:
            prompt = "(Cylinder) Visualize the VIS_FIELD of the cylinder to represent the 2D turbulent flow field of the object in the given image.";
            return cylinder(x, y, z, obj_center, float3(0u, 0u, lbm.get_Nz()), size);
        case 1:
            prompt = "(Sphere) Visualize the VIS_FIELD of the sphere to represent the 2D turbulent flow field of the object in the given image.";
            return sphere(x, y, z, obj_center, size);
        case 2:
            prompt = "(Ellipsoid) Visualize the VIS_FIELD of the ellipsoid to represent the 2D turbulent flow field of the object in the given image.";
            return ellipsoid(x, y, z, obj_center, float3(size, size, size/2));
        case 3:
            prompt = "(Cube) Visualize the VIS_FIELD of the cube to represent the 2D turbulent flow field of the object in the given image.";
            return cube(x, y, z, obj_center, size);
        case 4:
            prompt = "(Cuboid) Visualize the VIS_FIELD of the cuboid to represent the 2D turbulent flow field of the object in the given image.";
            return cuboid(x, y, z, obj_center, float3(size, size/2, size/2));
        case 5:
            prompt = "(Cone) Visualize the VIS_FIELD of the cone to represent the 2D turbulent flow field of the object in the given image.";
            return cone(x, y, z, obj_center, float3(0u, 0u, lbm.get_Nz()), size, size/2);
        case 6:
            prompt = "(Pipe) Visualize the VIS_FIELD of the pipe to represent the 2D turbulent flow field of the object in the given image.";
            return pipe(x, y, z, obj_center, float3(0u, 0u, lbm.get_Nz()), size);
        case 7:
            prompt = "(Cone Pipe) Visualize the VIS_FIELD of the cone pipe to represent the 2D turbulent flow field of the object in the given image.";
            return conepipe(x, y, z, obj_center, float3(0u, 0u, lbm.get_Nz()), size, size/2);
        case 8:
            prompt = "(Triangle) Visualize the VIS_FIELD of the triangle to represent the 2D turbulent flow field of the object in the given image.";
            return triangle(x, y, z, obj_center, float3(obj_center.x + size, obj_center.y, obj_center.z), float3(obj_center.x, obj_center.y + size, obj_center.z));
        case 9:
            prompt = "(Plane) Visualize the VIS_FIELD of the plane to represent the 2D turbulent flow field of the object in the given image.";
            return plane(x, y, z, obj_center, float3(0u, 1u, 0u));
        case 10:
            prompt = "(Torus X) Visualize the VIS_FIELD of the torus (X-axis) to represent the 2D turbulent flow field of the object in the given image.";
            return torus_x(x, y, z, obj_center, size/2, size);
        case 11:
            prompt = "(Torus Y) Visualize the VIS_FIELD of the torus (Y-axis) to represent the 2D turbulent flow field of the object in the given image.";
            return torus_y(x, y, z, obj_center, size/2, size);
        case 12:
            prompt = "(Torus Z) Visualize the VIS_FIELD of the torus (Z-axis) to represent the 2D turbulent flow field of the object in the given image.";
            return torus_z(x, y, z, obj_center, size/2, size);
        default:
            prompt = "Invalid ID.";
            return false;
    }
}

void run_simulation_2d(float si_u, int id, float size, bool is_train, bool render = true) {
    uint3 lbm_N = resolution(float3(1.0f, 2.0f, 0.01f), 1000u); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	// set the third dimension to 1u to make the simulation 2D
    lbm_N.z = 1u;
    const float si_length = 2.4f;
	// const float si_T = 10.0f;
	const float si_nu=1.48E-5f, si_rho=1.225f;
	const float lbm_length = 0.65f*(float)lbm_N.y; // length of the simulation box in lbm
	const float lbm_u = 0.055f; // velocity of the fluid in lbm
	units.set_m_kg_s(lbm_length, lbm_u, 1.0f, si_length, si_u, si_rho);        
	const float lbm_nu = units.nu(si_nu);
	// const ulong lbm_T = units.t(si_T);
	print_info("Re = "+std::to_string(to_uint(units.si_Re(si_length, si_u, si_nu))));
    // initialize the LBM object
	LBM lbm(lbm_N, lbm_nu);
    const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();
    size = size * Nx/4;
    float3 rotation = float3(0.0f, 0.0f, 0.0f);
    string prompt = id2object(id, rotation, lbm, size);
    // string prompt;

    // find the boundary location
    const uint Nx_1 = Nx - 1;
    const uint Ny_1 = Ny - 1;
    const uint Nz_1 = Nz - 1;
    parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		//if(z==0u) lbm.flags[n] = TYPE_S; // solid floor
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = lbm_u; // initialize y-velocity everywhere except in solid cells
		// if(id2object_2d(id, x, y, z, lbm, size, prompt)) lbm.flags[n] = TYPE_S;        
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_E; // all other simulation box boundaries are inflow/outflow
	}); // ####################################################################### run simulation, export images and data ##########################################################################
    // set the Visualize flag
    lbm.graphics.visualization_modes = VIS_FIELD|VIS_FLAG_LATTICE;
    lbm.graphics.slice_mode = 3;
    if (render){

        float dt = units.si_t(1ul); // 1 LBM步对应的真实时间，单位是秒
        print_info("dt = " + std::to_string(dt));
        float slow_down = 100.0f; 
        int frequency = (int)((1/24.0)/dt/slow_down); // 24fps
        // int frequency = 1;
        print_info("frequency = " + std::to_string(frequency));
        create_dataset_2d(frequency, is_train, lbm, prompt);
    }else{
        
        lbm.run();
    }
}
