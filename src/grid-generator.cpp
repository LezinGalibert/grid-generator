#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <map>
#include <cstdlib>

float layer_property_normal_distribution(float mean, float std_deviation) {

	std::default_random_engine generator;
	std::normal_distribution<float> distribution(mean, std_deviation);

	return distribution(generator);
}

float compute_permeability(float porosity, float Amin, float Amax, float b) {

	  std::default_random_engine generator;
	  std::uniform_real_distribution<double> distribution(0.0,1.0);

	  float alpha = distribution(generator);

	  return pow(porosity, b)*(alpha*Amin + (1 - alpha)*Amax);
}

float compute_normal_sat(float So, float Sf, float S) {

	return (S - So)/(Sf - So);
}

float compute_inv_normal_sat(float So, float Sf, float Sn) {

	return Sn*(Sf - So) + So;
}

//kr =  a*exp(-lambda*S) + b
// with kr(So) = kri
// with kr(Sf) = 0
float compute_relperm(float Sn, float N, float kri = 1) {

	return kri*pow(Sn, N);
}

//J = aS^(-b)
float compute_power_J(float Sn, float a, float b) {

	return a*pow(Sn, -b);
}

float compute_inv_power_J(float J, float a, float b) {

	return pow(J/a, -1/b);
}

//J = 1/S - 1/(1-S)
float compute_sym_J(float Se) {

	return 1/Se - 1/(1-Se);
}

float compute_inv_sym_J(float J) {

	if (J != 0) {
		return (2 + J - pow(J*J + 4, 0.5))/(2*J);
	}
	else {
		return 0.5;
	}
}

float compute_weighted_average(float k1, float k2, float ratio, int exponent) {

	return pow(ratio*pow(k1, exponent) + (1 - ratio)*pow(k2, exponent), 1/exponent);
}

std::string build_saturation_table(std::vector<float> water_saturation,
		                           std::vector<float> kr_water,
		                           std::vector<float> kr_oil,
		                           std::vector<float> Jfunction) {

	size_t size {water_saturation.size()};
	std::string table_row {""};

	if ((kr_water.size()  != size) |
		(kr_oil.size()    != size) |
		(Jfunction.size() != size)) {

		std::cerr << "ERROR: Columns have different sizes." << std::endl;
	}

	for (size_t row {0}; row < size; ++row) {

		table_row +=        std::to_string(water_saturation[row]) +
				     "\t" + std::to_string(kr_water[row])         +
					 "\t" + std::to_string(kr_oil[row])           +
					 "\t" + std::to_string(Jfunction[row])        + "\n";
	}

	return table_row;
}

template <class Container>
void split(const std::string& str, Container& cont,
           std::string delim)
{
    std::size_t current, previous = 0;
    current = str.find(delim);
    while (current != std::string::npos) {
        cont.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    cont.push_back(str.substr(previous, current - previous));
}

std::map<std::string, std::vector<std::string>> parse_file(std::string filename){

	std::ifstream infile;
	std::map<std::string, std::vector<std::string>> map;
	infile.open(filename);

	std::string line;
	while (std::getline(infile, line))
	{
	    std::vector<std::string> split_line;

	    split(line, split_line, "=");
	    split(split_line[1], map[split_line[0]], " ");
	}

	return map;
}

std::vector<int> generate_layer_map(int pattern_layer_number, float thickness_ratio, int depth, int rock_number = 2) {

	// 2 rocks, 3 layers, ratio of 0.75, depth of 16

	//    ****************
	//    ****************
	//    ****************
	//    ****************
	//    ++++++++++++++++
	//    ++++++++++++++++
	//    ****************
	//    ****************
	//    ****************
	//    ****************
	//    ++++++++++++++++
	//    ++++++++++++++++
	//    ****************
	//    ****************
	//    ****************
	//    ++++++++++++++++

	int pattern_thickness {depth/pattern_layer_number};
	std::vector<int> rock_thickness {pattern_thickness*thickness_ratio, pattern_thickness*(1 - thickness_ratio)};
	int thickness_leftover {depth - pattern_layer_number*(rock_thickness[0] + rock_thickness[1])};

	int current_depth {0};
	int current_rock {1};

	int layer_limit {};

	std::vector<int> layer_map {};

	while (current_depth < depth ) {

		layer_limit = (current_depth < depth - thickness_leftover) ? rock_thickness[current_rock - 1] : thickness_leftover;

		for (int d {0}; d < layer_limit; ++d){

			layer_map.push_back(current_rock);
		}
		current_depth += rock_thickness[current_rock - 1];

		if (thickness_leftover > 0) {
			layer_map.push_back(current_rock);
			current_depth++;
			thickness_leftover--;
		}

		current_rock %= rock_number;
		current_rock += 1;
	}

	return layer_map;
}

int main(int argc, char **argv) {

	std::string grid_file_name {argv[1]};
	std::map<std::string, std::vector<std::string>> param_map;
	std::vector<int> layer_mask;
	param_map = parse_file("/home/lezin/OPM/dataset/petrophysics/" + grid_file_name);

	int rock_number {stoi(param_map["rock_number"][0])};
	int layer_number {stoi(param_map["layer_number"][0])};
	float thickness_ratio {stof(param_map["thickness_ratio"][0])};

	std::string path_to_folder {"/home/lezin/OPM/dataset/" + param_map["case_name"][0] + '/'};

	float x_min {stof(param_map["x_min"][0])};
	float x_max {stof(param_map["x_max"][0])};
	float y_min {stof(param_map["y_min"][0])};
	float y_max {stof(param_map["y_max"][0])};
	float z_min {stof(param_map["z_min"][0])};
	float z_max {stof(param_map["z_max"][0])};

	int nx {stoi(param_map["nx"][0])};
	int ny {stoi(param_map["ny"][0])};
	int nz {stoi(param_map["nz"][0])};

	float sigma {stof(param_map["sigma"][0])};

	float std_deviation {stof(param_map["std_deviation"][0])};

	float So {stof(param_map["So"][0])};
	float Sf {stof(param_map["Sf"][0])};

	std::ofstream param;
	std::ofstream sat_pd;
	std::ofstream rock_list;

	float stepx {0};
	float stepy {0};
	float stepz {0};

	std::vector<float> krio {stof(param_map["rock_krio"][0]), stof(param_map["rock_krio"][1])};
	std::vector<float> Nw {stof(param_map["Nw"][0]), stof(param_map["Nw"][1])};
	std::vector<float> No {stof(param_map["No"][0]), stof(param_map["No"][1])};
	std::vector<float> a {stof(param_map["rock_a"][0]), stof(param_map["rock_a"][1])};
	std::vector<float> b {stof(param_map["rock_b"][0]), stof(param_map["rock_b"][1])};

	std::vector<float> mean_poro {stof(param_map["mean_porosities"][0]), stof(param_map["mean_porosities"][1])};
	std::vector<float> A_min {stof(param_map["A_min"][0]), stof(param_map["A_min"][1])};
	std::vector<float> A_max {stof(param_map["A_max"][0]), stof(param_map["A_max"][1])};
	std::vector<float> corr {stof(param_map["corr"][0]), stof(param_map["corr"][1])};

	std::string flow_direction {param_map["flow_direction"][0]};
	bool parallel_flow {flow_direction != "2"};

	int num_sats {stoi(param_map["num_sats"][0])};
	float min_sat {stof(param_map["min_sat"][0])};
	float max_sat {stof(param_map["max_sat"][0])};

	int num_pdrops {stoi(param_map["num_pdrops"][0])};
	float min_pdrop {stof(param_map["min_pdrop"][0])};
	float max_pdrop {stof(param_map["max_pdrop"][0])};

	std::string specgrid {"SPECGRID\n"};
	std::string coord {"COORD\n"};
	std::string zcorn {"ZCORN\n"};
	std::string satnum {"SATNUM\n"};
	std::string poro {"PORO\n"};
	std::string permx {"PERMX\n"};
	std::string permy {"PERMY\n"};
	std::string permz {"PERMZ\n"};
	std::string tmp {""};

	std::default_random_engine generator;

	float Sw1, Sw2, J1, J2, krw1, krw2, kro1, kro2, krwe, kroe;

	param.open(path_to_folder + "parameters.param");
	param << "fileformat=eclipse\nfilename="
		  << path_to_folder
		  << "eclipse_grid.grdecl\nrock_list="
		  << path_to_folder
		  << "rocklist.txt\nboundary_condition_type=" + param_map["boundary_condition_type"][0]
		  << "\nflow_direction=" + flow_direction
		  << "\nsigma=" + std::to_string(sigma);
	param.close();


	sat_pd.open(path_to_folder + "sat_pdrop.param");
	sat_pd << num_sats << std::endl;
	for (float sat {min_sat}; sat < max_sat; sat += (max_sat - min_sat)/num_sats){

		sat_pd << " " << std::fixed << std::to_string(sat) << " " << std::to_string(num_pdrops);

		for (float press {min_pdrop}; press <= max_pdrop; press *= pow(max_pdrop/min_pdrop, 1.0/num_pdrops)){
			sat_pd << " " << std::scientific << std::to_string(press);

			if (num_pdrops == 1) break;
		}
		sat_pd << std::endl;
	}
	sat_pd.close();


	layer_mask = generate_layer_map(layer_number, thickness_ratio, nz, rock_number);

	rock_list.open(path_to_folder + "rocklist.txt");
	rock_list << std::to_string(rock_number);
	for (int r {1}; r <= rock_number; ++r) {
		rock_list << "\nrocktypes/rock" << std::to_string(r) << ".txt";
	}
	rock_list.close();


	specgrid += " " + std::to_string(nx);
	specgrid += " " + std::to_string(ny);
	specgrid += " " + std::to_string(nz);
	specgrid += " 1 F /";


	stepx = (x_max - x_min)/nx;
	stepy = (y_max - y_min)/ny;
	stepz = (z_max - z_min)/nz;

	for (int j {0}; j <= ny; ++j) {
		for (int i {0}; i <= nx; ++i) {

			tmp = " " + std::to_string(i*stepx + x_min);
			tmp += " " + std::to_string(j*stepy + y_min);
			coord += tmp + " " + std::to_string(z_min) + tmp + " " + std::to_string(z_max) + "\n";
		};
	};

	coord += "/";


	zcorn += " " + std::to_string(4*nx*ny) + "*" + std::to_string(z_min) + "\n";
	tmp = "";

	for (int k {1}; k < nz; ++k) {
		tmp = " " + std::to_string(4*nx*ny) + "*" + std::to_string(z_min + k*stepz) + "\n";
		zcorn += tmp + tmp;
	};

	zcorn += " " + std::to_string(4*nx*ny) + "*" + std::to_string(z_max) + "\n";
	zcorn += "/";


	for (int r {0}; r < rock_number; ++r){

		float Sn, Se;

		std::vector<float> *water_sats = new std::vector<float>;
		std::vector<float> *kr_water = new std::vector<float>;
		std::vector<float> *kr_oil = new std::vector<float>;
		std::vector<float> *Jfunction = new std::vector<float>;

		for (float S {So}; S <= Sf; S+= So) {

			Sn = compute_normal_sat(0, 1, S);
			Se = compute_normal_sat(0, 1, S);
			water_sats->push_back(S);
			kr_water->push_back(compute_relperm(Sn, Nw[r]));
			kr_oil->push_back(compute_relperm(1-Sn, No[r], krio[r]));
			Jfunction->push_back(compute_sym_J(Se));
		};


		std::string *saturations = new std::string {"-- Sw Krw Kro J\n"};

		try
		{
			*saturations += build_saturation_table(*water_sats, *kr_water, *kr_oil, *Jfunction);
		}
		catch (int e)
		{
			std::cout << "An exception occurred. Exception Nr. " << e << '\n';
		};

		delete Jfunction;
		delete kr_oil;
		delete kr_water;
		delete water_sats;

        //std::cout << *saturations << "\n" << std::endl;


		std::ofstream *satfile = new std::ofstream;

		satfile->open(path_to_folder + "rocktypes/rock" + std::to_string(r + 1) + ".txt");
		*satfile << *saturations;
		satfile->close();

		delete satfile;
		delete saturations;

	}

	float sampled_poro, computed_permx;

	for (std::vector<int>::iterator rock {layer_mask.begin()}; rock != layer_mask.end(); rock++){

		std::normal_distribution<float> distribution(mean_poro[*rock - 1], mean_poro[*rock - 1]*std_deviation);

		satnum += ' ' + std::to_string(nx*ny) + '*' + std::to_string(*rock) + "\n";

		for (int j {0}; j < ny; ++j){
			for (int i {0}; i < nx; ++i){

				sampled_poro = distribution(generator);
				poro += " " + std::to_string(sampled_poro);

				computed_permx = compute_permeability(sampled_poro, A_min[*rock - 1], A_max[*rock - 1], corr[*rock - 1]);
				permx += " " + std::to_string(computed_permx);
			}
			poro += "\n";
			permx += "\n";
		}
	}


	satnum += "/";
	poro += "/";
	permx += "/";
	permy += "/";
	permz += "/";

//	std::cout << specgrid << "\n "<< std::endl;
//	std::cout << coord << "\n " << std::endl;
//	std::cout << zcorn << "\n " << std::endl;
//	std::cout << satnum << "\n " << std::endl;
//	std::cout << poro << "\n " << std::endl;
//	std::cout << permx << "\n " << std::endl;

	std::ofstream outfile;

	outfile.open(path_to_folder + "eclipse_grid.grdecl");
	outfile << specgrid << "\n "<< std::endl;
	outfile << coord << "\n " << std::endl;
	outfile << zcorn << "\n " << std::endl;
	outfile << satnum << "\n " << std::endl;
	outfile << poro << "\n " << std::endl;
	outfile << permx << "\n " << std::endl;

	outfile.close();

	std::vector<float> capillary_pressures;
	std::vector<float> eq_water_sats;
	std::vector<float> eq_water_relperm;
	std::vector<float> eq_oil_relperm;

	for (float cp {-200000000}; cp <= 500000000; cp+= 1000){

		J1 = cp/sigma*pow(A_min[0]*9.869233E-16/mean_poro[0], 0.5);
		J2 = cp/sigma*pow(A_min[1]*9.869233E-16/mean_poro[1], 0.5);

		Sw1 = compute_inv_sym_J(J1);
		Sw2 = compute_inv_sym_J(J2);

		krw1 = compute_relperm(Sw1, Nw[0]);
		krw2 = compute_relperm(Sw2, Nw[1]);

		kro1 = compute_relperm(1 - Sw1, No[0]);
		kro2 = compute_relperm(1 - Sw2, No[1]);

//		std::cout << "\n**********" << "\nCP=" << std::to_string(cp)
//				                    << "\nJ1=" << std::to_string(J1)
//				                    << "\nJ2=" << std::to_string(J2)
//									<< "\nSw1=" << std::to_string(Sw1)
//									<< "\nSw2=" << std::to_string(Sw2)
//		                            << "\nSeq=" << std::to_string(thickness_ratio*Sw1 + (1 - thickness_ratio)*Sw2)
//									<< "\nkrw1=" << std::to_string(krw1)
//									<< "\nkrw2=" << std::to_string(krw2)
//									<< "\nkro1=" << std::to_string(kro1)
//									<< "\nkro2=" << std::to_string(kro2);

		if (parallel_flow) {

			krwe  = compute_weighted_average(krw1*A_min[0], krw2*A_min[1], thickness_ratio, 1);
			krwe /= compute_weighted_average(     A_min[0],      A_min[1], thickness_ratio, 1);

			kroe  = compute_weighted_average(kro1*A_min[0], kro2*A_min[1], thickness_ratio, 1);
			kroe /= compute_weighted_average(     A_min[0],      A_min[1], thickness_ratio, 1);

		} else {

			krwe  = compute_weighted_average(krw1*A_min[0], krw2*A_min[1], thickness_ratio, -1);
			krwe /= compute_weighted_average(     A_min[0],      A_min[1], thickness_ratio, -1);

			kroe  = compute_weighted_average(kro1*A_min[0], kro2*A_min[1], thickness_ratio, -1);
			kroe /= compute_weighted_average(     A_min[0],      A_min[1], thickness_ratio, -1);

		}

		capillary_pressures.push_back(cp);
		eq_water_sats.push_back(thickness_ratio*Sw1 + (1 - thickness_ratio)*Sw2);
		eq_water_relperm.push_back(krwe);
		eq_oil_relperm.push_back(kroe);
	};

	std::string *eq_saturations = new std::string {"-- Sw Krw Kro Pc\n"};

	try
	{
		*eq_saturations += build_saturation_table(eq_water_sats, eq_water_relperm, eq_oil_relperm, capillary_pressures);
	}
	catch (int e)
	{
		std::cout << "An exception occurred. Exception Nr. " << e << '\n';
	};

//	std::cout << *eq_saturations;

	std::ofstream outfile_eq_sats;

	outfile_eq_sats.open(path_to_folder + "analytical_relperms");
	outfile_eq_sats << *eq_saturations;
	outfile_eq_sats.close();

	delete eq_saturations;

	std::ofstream outfile_capillary_number;
	outfile_capillary_number.open(path_to_folder + "capillary_number.params");
	outfile_capillary_number << "Flow direction\t" << param_map["flow_direction"][0] << std::endl;
	outfile_capillary_number << "Kt\t";

	if (parallel_flow) {
		outfile_capillary_number << compute_weighted_average(A_min[0], A_min[1], thickness_ratio,  1);
	}
	else {
		outfile_capillary_number << compute_weighted_average(A_min[0], A_min[1], thickness_ratio, -1);
	}
	outfile_capillary_number << std::endl;
	outfile_capillary_number << "Sigma\t" << sigma << std::endl;
	outfile_capillary_number << "L\t";

	if (flow_direction == "0") {
		outfile_capillary_number << std::to_string(x_max - x_min);
	}
	else if (flow_direction == "1") {
		outfile_capillary_number << std::to_string(y_max - y_min);
	}
	else if (flow_direction == "2") {
		outfile_capillary_number << std::to_string(z_max - z_min);
	}
	outfile_capillary_number << std::endl;
	outfile_capillary_number << "Thickness ratio\t" << std::to_string(thickness_ratio);

	outfile_capillary_number << std::endl;


	return 0;
}
