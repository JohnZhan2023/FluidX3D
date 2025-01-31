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
            shape_name = "cylinder";
            break;
        case 1:
            shape_name = "sphere";
        
            break;
        case 2:
            shape_name = "ellipsoid";
            break;
        case 3:
            shape_name = "cube";
            break;
        case 4:
            shape_name = "cone";
            break;
        case 5:
            shape_name = "pipe";
            break;
        case 6:
            shape_name = "cone pipe";
            break;
        case 7:
            shape_name = "cuboid";
            break;
        case 8:
            shape_name = "triangle";
            break;
        case 9:
            shape_name = "torus (X axis)";
            break;
        case 10:
            shape_name = "torus (Y axis)";
            break;
        case 11:
            shape_name = "torus (Z axis)";
            break;
        case 12:
            shape_name = "dog";
            lbm.voxelize_stl(get_exe_path()+"../stl/Dog_Singlecolor.stl", obj_center, rotation, size);
            break;
        case 13:
            shape_name = "shark";
            lbm.voxelize_stl(get_exe_path()+"../stl/jeff_the_land_shark.stl", obj_center, rotation, size);
            break;
        case 14:
            shape_name = "landing gear";
            lbm.voxelize_stl(get_exe_path()+"../stl/crm-hl_reference_ldg.stl", obj_center, rotation, size);
            break;
        case 15:
            shape_name = "cat ring";
            lbm.voxelize_stl(get_exe_path()+"../stl/cat ring.stl", obj_center, rotation, size);
            break;
        case 16:
            shape_name = "flexible cat";
            lbm.voxelize_stl(get_exe_path()+"../stl/CatFlexi.stl", obj_center, rotation, size);
            break;
        case 17:
            shape_name = "gong-gi";
            lbm.voxelize_stl(get_exe_path()+"../stl/gong-gi.stl", obj_center, rotation, size);
            break;
        case 18:
            shape_name = "octopus";
            lbm.voxelize_stl(get_exe_path()+"../stl/Octopus_spiral_v6.stl", obj_center, rotation, size);
            break;
        case 19:
            shape_name = "planter";
            lbm.voxelize_stl(get_exe_path()+"../stl/radiant_planter.stl", obj_center, rotation, size);
            break;
        case 20:
            shape_name = "plate";
            lbm.voxelize_stl(get_exe_path()+"../stl/radiant_plate.stl", obj_center, rotation, size);
            break;
        case 21:
            shape_name = "shark v2.0";
            lbm.voxelize_stl(get_exe_path()+"../stl/Shark_v2.0_A.stl", obj_center, rotation, size);
            break;
        case 22:
            shape_name = "car wheel SKODA";
            lbm.voxelize_stl(get_exe_path()+"../stl/SKODA_ENYAQ_FL_COUPE_S24_WHEEL_C5W_R21_20241220.stl", obj_center, rotation, size);
            break;
        case 23:
            shape_name = "Zorro character";
            lbm.voxelize_stl(get_exe_path()+"../stl/zorro_chibi2fix.stl", obj_center, rotation, size);
            break;
        case 24:
            shape_name = "Cow";
            lbm.voxelize_stl(get_exe_path()+"../stl/Cow_t.stl", obj_center, rotation, size);
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
    string original_path = "C:\\Users\\86177\\Desktop\\workspace\\EXP\\data_generation\\";

    const uint star_T = 1000u; // number of LBM time steps to simulate
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
    float3(((float)Nz/Nx) * Nx, 0.0f * Ny, 0.0f * Nz),  // 相机位置（x, y, z）
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
				float3(((float)Nz/Nx) * Nx, 0.0f * Ny, 0.0f * Nz), // 相机位置（x, y, z）
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
        int frequency = (int)((1/24.0)/dt);
        // int frequency = 1;
        print_info("frequency = " + std::to_string(frequency));
        create_dataset(frequency, is_train, lbm, prompt);
    }else{
        lbm.run();
    }

}

