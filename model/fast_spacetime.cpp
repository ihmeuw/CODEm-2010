#include <iterator>
#include <iostream>
#include <cstdlib>
#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>

using namespace std ;

#define AGE_EPS 0.0001

class AgeIndex {
            public:
				AgeIndex() ;
                AgeIndex(double start_age,  double end_age) ;  // constructor
				AgeIndex(AgeIndex &input) ; //copy constructor
				void initialize() ; // build age_index
				void deallocate() ; // free memory
                ~AgeIndex() {}    // destructor

				int numGroups() const; // return number of age groups
				int getIndex(double age) const; // retrieve index given an age

            protected:
				int num_age_groups ; // number of age groups
				vector<int> age_index ;     // Vector to store age and corresponding indicies
				double start_age ;
				double end_age ;
        };

        //Constructor
		AgeIndex::AgeIndex() {} ;

        AgeIndex::AgeIndex(double start_age, double end_age) : start_age(start_age),end_age(end_age) {};
		void AgeIndex::initialize() {

			// To avoid an O(n) lookup for every age group's index
			// I made the total array much longer, so the index of each age group
			// Is located in age_index[age_group*100]

			// Age group only needs to be as large as end_age*100
			// Add 0.5 for floating point arithmetic
			// Forcing to be int will round down

            cout << "end_age: " << end_age  << endl;
            cout << "int(end_age*100 + 0.5): " << int(end_age*100 + 0.5) << endl;
			age_index.resize(int(end_age*100 + 1.5)) ;


			// Enumerate age groups for weighting and iterating
			// This method assumes that age groups are taken from the set:
			// 0,0.01,0.1,1,5,10,15,...,75,80
			double age_difference= end_age - start_age ;


			if (age_difference >= 1 || end_age == 1) {
				num_age_groups = (int)floor(age_difference/5 + 0.5) ;
				if (start_age == 0)
					num_age_groups += 4 ;
				else if (start_age == 0.01)
					num_age_groups += 3 ;
				else if (start_age == 0.1)
					num_age_groups += 2 ;
				else
					num_age_groups += 1 ;
			}
			// 0 to 0.1
			else if (fabs(age_difference- 0.1) < AGE_EPS)
				num_age_groups = 3 ;
			// 0 to 0.01
			else if (fabs(age_difference- 0.01) < AGE_EPS)
				num_age_groups = 2 ;
			// 0.01 to 0.1
			else if (fabs(age_difference- 0.09) < AGE_EPS)
				num_age_groups = 2 ;
			// Same age, not sure if this ever happens
			else if (fabs(age_difference- 0) < AGE_EPS)
				num_age_groups = 1 ;


			// fill in indicies
			double current_age = start_age ;
			bool passed_1 = (start_age > 1) ;
			for (int current_group = 0; current_group <num_age_groups; current_group++) {
				if (!passed_1) {
					if (fabs(current_age-0) < AGE_EPS) {
						age_index[0] = current_group ;
						current_age = 0.01 ;
					}
					else if (fabs(current_age-0.01) < AGE_EPS) {
						age_index[1] = current_group ;
						current_age = 0.1 ;
					}
					else if (fabs(current_age-0.1) < AGE_EPS) {
						age_index[10] = current_group ;
						current_age = 1;
					}
					else if (fabs(current_age-1) < AGE_EPS) {
						age_index[100] = current_group ;
						current_age = 5;
						passed_1 = true ;
					}
				}
				else {
					age_index[int(current_age*100 + 0.5)] = current_group;
					current_age = current_age + 5 ;
				}
			}
		}

		void AgeIndex::deallocate() {

		}


		//For a given age, return the index
		int AgeIndex::getIndex(double age) const {
			//char display_index[100] ;
			//sprintf(display_index,"Inside getIndex: %f", end_age) ;
			//SF_display(display_index) ;
			//SF_display("\n") ;
			return age_index[int(age*100 + 0.5)];
		}

		//Return number of age groups
		int AgeIndex::numGroups() const{
			//char display_groups[100] ;
			//sprintf(display_groups,"Number of Age Groups: %d",num_age_groups) ;
			//SF_display(display_groups) ;
			return num_age_groups;
		}

class AgeWeights {
            public:

                AgeWeights (const AgeIndex &age_index, double omega): age_index(age_index), omega(omega) {} ;  // constructor
                ~AgeWeights() {}    // destructor

				double GetWeight(double age_1, double age_2) const ;

				void initialize() ; // build age_index
				void deallocate() ; // free memory



            protected:

				vector<double> age_weights; // 1Dvector to store age weights, representing an 2d matrix
				// num_age_groups by num_age_groups
				int num_age_groups; //stored here for convenience, from age_index
				const AgeIndex& age_index ;
				double omega;
        };

        //Constructor


		void AgeWeights::initialize () {
			num_age_groups = age_index.numGroups() ;
			// Create age weights array

			age_weights.resize(  (int) ( (num_age_groups*(num_age_groups+1))/2 +0.5), 0.0) ;

			for (int i = 0; i < num_age_groups; ++i) {
				for (int j = 0; j < i+1; j++) {
					age_weights[((i+1)*i)/2+j] = 1/(exp(omega*fabs(double(i-j)))) ;
				}
			}

		}

		//For two ages, return the weight
		double AgeWeights::GetWeight(double age_1, double age_2) const {
			if (age_1 <= age_2) {
				return age_weights[((age_index.getIndex(age_2)+1)*age_index.getIndex(age_2))/2+age_index.getIndex(age_1)];
			}
			else {
				return age_weights[((age_index.getIndex(age_1)+1)*age_index.getIndex(age_1))/2+age_index.getIndex(age_2)];
			}
		}

		void AgeWeights::deallocate() {		  }

class YearWeights {
			public:
				YearWeights(int start_year, int end_year): start_year(start_year), end_year(end_year) {} ; //constructor
				~YearWeights() {} // destructor
				double GetWeight(double year_1, double year_2, double lambda) const ;
				void initialize () ;
			protected:

				int start_year;
				int end_year;
} ;

		void YearWeights::initialize () {

		}

		double YearWeights::GetWeight(double year_1, double year_2, double lambda) const {
			double numerator = fabs(year_1 - year_2) ;
			double denominator = max(fabs(year_1-start_year), fabs(year_1-end_year)) + 1.0 ;
			double ratio = numerator / denominator ;
			ratio  = pow(ratio, lambda) ;

			return pow(1-ratio, 3) ;
		}


// csv parsing code from http://stackoverflow.com/a/1120224/1748679
class CSVRow
{
    public:
        string const& operator[](size_t index) const
        {
            return m_data[index];
        }
        size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(istream& str)
        {
            string         line;
            getline(str,line);

            stringstream   lineStream(line);
            string         cell;

            m_data.clear();
            while(getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        vector<string>    m_data;
};

istream& operator>>(istream& str,CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

class CSVIterator
{
    public:
        typedef input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }
        CSVIterator()                   :m_str(NULL) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { (*m_str) >> m_row;m_str = m_str->good()?m_str:NULL;}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        istream*       m_str;
        CSVRow              m_row;
};


	struct weight {
		int relatedness ;
		double wt;
	};

	struct country_age {
		string iso3 ;
		double age ;

	};
bool operator==(const country_age& l, const country_age& r)
		{
			return ( (l.iso3 == r.iso3) && (fabs(l.age - r.age) < AGE_EPS) ) ;
		};
int main(int argc, char* argv[])
{
	cout << "BEFORE SET THREADS" << endl;
	omp_set_num_threads(12);
	cout <<"START NOW" << endl;
	int num_submodels = atoi(argv[1]);
	double start_age = atof(argv[2]) ;
	double end_age = atof(argv[3]) ;
	int start_year = atoi(argv[4]) ;
	int end_year = atoi(argv[5]) ;
	double lambda = atof(argv[6]) ;
	double lambda_no_data = atof(argv[7]) ;
	double zeta = atof(argv[8]) ;
	double zeta_no_data = atof(argv[9]) ;
	double omega = atof(argv[10]) ;
	string test = argv[11];
	string inname = argv[12];
	string outname = argv[13] ;


    // Make threads race to print for fun and to prove you are parallel
	#pragma omp parallel for
	for(int i = 0; i < argc; i++)
      cout << "argv[" << i << "] = " << argv[i] << endl;

    // Get number of observations
    ifstream       file2(inname.c_str());

    CSVIterator get_rows(file2);
    int rows = 0;
    ++get_rows; // header row
    for(get_rows; get_rows != CSVIterator(); ++get_rows,++rows) {
    }
    int num_obs = rows;
    cout << "rows: " << num_obs << endl;


    ifstream       file(inname.c_str());
	CSVIterator loop(file) ;
	cout << "Woah there, I read opened the file" << endl;
	cout << loop->size() ;

	// Column where variable is stored
	int iso3_loc;
	int age_loc;
	int year_loc;
	int region_loc;
	int super_region_loc;
	int national_loc;
	int has_data_loc;
	int ln_rate_loc ;
	int predictme_loc;
	vector<int> linear_locs (num_submodels);
	vector<int> residual_locs (num_submodels);


	#pragma omp parallel for
	for(int i = 0; i < loop->size() ; ++i) {
		if ((*loop)[i] == "iso3")
			iso3_loc = i;

		if ((*loop)[i] == "age")
			age_loc = i;

		if ((*loop)[i] == "year")
			year_loc = i;

		if ((*loop)[i] == "region")
			region_loc = i;

		if ((*loop)[i] == "super_region")
			super_region_loc = i;

		if ((*loop)[i] == "national")
			national_loc = i;

		if ((*loop)[i] == "has_data")
			has_data_loc = i;
		if ((*loop)[i] == "ln_rate")
			ln_rate_loc = i;

		if ((*loop)[i] == "predictme_"+test)
			predictme_loc = i;

		if ((*loop)[i].substr(0,6) == "linear"){

				// (*loop)[i].find_first_of(_) gets index of underscore
				// i.e. the "_" in linear_123
				// (*loop)[i].substr( non_numeric+1) returns the number at the end of the string
				// i.e. 123 in linear_123
				// atoi (.....c_str) converts this to an int, and 1 is subtracted to index from 0 to
				// number_submodels -1
				int underscore_index = (*loop)[i].find_first_of("_") ;
				string submodel_as_string = (*loop)[i].substr(underscore_index + 1) ;
				int submodel_as_int = atoi(submodel_as_string.c_str()) ;
				linear_locs[submodel_as_int - 1] = i ;

		}

		if ((*loop)[i].substr(0,8) == "residual"){
				// (*loop)[i].find_last_not_of(01....) gets index of last non-numeric char
				// i.e. the "_" in linear_123
				// (*loop)[i].substr( non_numeric+1) returns the number at the end of the string
				// i.e. 123 in linear_123
				// atoi (.....c_str) converts this to an int, and 1 is subtracted to index from 0 to
				// number_submodels -1
				int underscore_index = (*loop)[i].find_first_of("_") ;
				string submodel_as_string = (*loop)[i].substr(underscore_index + 1) ;
				int submodel_as_int = atoi(submodel_as_string.c_str()) ;
				residual_locs[submodel_as_int - 1 ] = i ;

		}
	}

	// Move past the header row
	++loop;

	// Create way to store data
	vector<string> iso3(num_obs);
	vector<double> age(num_obs);
	vector<int> year(num_obs);
	vector<int> region(num_obs);
	vector<int> super_region(num_obs);
	vector<int> has_data(num_obs); //int to compute sum
	vector<int> predictme(num_obs);
	vector<int> national(num_obs);
	vector<int> has_data_rows; //has data
	vector<double> ln_rate(num_obs);
	vector<country_age> prediction_pairs ; //vector of country ages that get predicted
	vector<country_age> has_data_pairs; // vector of pairs with some data

	vector<double> lin_data(num_obs*num_submodels);
	vector<double> res_data(num_obs*num_submodels);


	int row_num = 0;
	int total_data = 0;
	int has_data_total = 0;
	cout << "before read" << endl;
	//openmp struggles with parallelizing this, so I didnt
	for(loop; loop != CSVIterator(); ++loop,++row_num) {

		iso3[row_num] = (*loop)[iso3_loc] ;
		age[row_num] = atof((*loop)[age_loc].c_str()) ;
		year[row_num] = atoi((*loop)[year_loc].c_str()) ;
		region[row_num] = atoi((*loop)[region_loc].c_str()) ;
		super_region[row_num] = atoi((*loop)[super_region_loc].c_str()) ;
		national[row_num] = atoi((*loop)[super_region_loc].c_str());
		has_data[row_num] = atoi((*loop)[has_data_loc].c_str()) ;
		if (has_data[row_num]) ++has_data_total;
		ln_rate[row_num] = atof((*loop)[ln_rate_loc].c_str()) ;
		predictme[row_num] = atoi((*loop)[predictme_loc].c_str()) ;
		if ( has_data[row_num] || predictme[row_num]) {
			country_age current_pair;
			current_pair.iso3 = iso3[row_num];
			current_pair.age = age[row_num];
			if (std::find(prediction_pairs.begin(), prediction_pairs.end(), current_pair)==prediction_pairs.end()) {
				prediction_pairs.push_back(current_pair);
			}
			if (has_data[row_num]) {
				if (std::find(has_data_pairs.begin(), has_data_pairs.end(), current_pair)==has_data_pairs.end()) {
					has_data_pairs.push_back(current_pair);
				}
				has_data_rows.push_back(row_num) ;
			}
		}

		for (int i = 0;i < num_submodels; ++i) {
			lin_data[row_num*num_submodels+i] =  atof((*loop)[linear_locs[i]].c_str()) ;
			res_data[row_num*num_submodels+i] =  atof((*loop)[residual_locs[i]].c_str()) ;
		}
	}
	cout << "Rows in CSV: " << row_num << endl;
	AgeIndex age_index(start_age, end_age);
	age_index.initialize() ;
	const AgeIndex& my_const_ref = age_index;
	AgeWeights age_weights(my_const_ref, omega) ;
	age_weights.initialize();

	YearWeights year_weights(start_year, end_year) ;
	year_weights.initialize();

	//Lets build the weight matrix

	// get the predictions we actually need to make
	vector<int> prediction_rows ;
	// 1 if there is data anywhere in country-age pair
	vector<int> has_data_pairwise(num_obs,0) ;

	for (int i = 0; i < num_obs; ++i) {
		country_age current_pair;
		current_pair.iso3 = iso3[i];
		current_pair.age = age[i];
		if (std::find(prediction_pairs.begin(), prediction_pairs.end(), current_pair)!=prediction_pairs.end()) {
				prediction_rows.push_back(i);
			}
		if (std::find(has_data_pairs.begin(), has_data_pairs.end(), current_pair)!=has_data_pairs.end()) {
				has_data_pairwise[i] = 1;
			}
	}




	vector<double> pred_residual(num_submodels*num_obs, 0.0) ;
	vector<double> simple_spacetime_ln_rate(num_obs, 0.0 ) ;

	cout << "Prediction Rows Size: "<< prediction_rows.size() << endl;
	cout << "Has data size: "<< has_data_rows.size() << endl;
	vector<double> residual_weighted_average_totals(prediction_rows.size(), 0.0)  ;
	vector<double> simple_weighted_average_totals(prediction_rows.size(), 0.0)  ;

	#pragma omp parallel for
	for (int i = 0; i < prediction_rows.size(); ++i) {
	//for (int i = 0; i < num_obs; ++i) {
		vector<weight> w_mat(has_data_rows.size()) ;

		//cout << i << endl;

		int row_num = prediction_rows[i] ;
		//int row_num = i;
		double age_i = age[row_num];
		int year_i = year[row_num];
		string country_i = iso3[row_num];
		int region_i = region[row_num];
		int super_region_i = super_region[row_num];
		double this_lambda;
		double this_zeta;
		if (  has_data_pairwise[row_num]) {
			this_lambda = lambda ;
			this_zeta = zeta;
		}
		else {
			this_lambda = lambda_no_data ;
			this_zeta = zeta_no_data ;
		}

		vector<double> row_sums(4, 0.0) ;
		double subnational_sum;
		double total_sum;


		for (int j=0; j < has_data_rows.size(); ++j) {
			int row_j = has_data_rows[j];
			double age_j = age[row_j];
			int year_j = year[row_j];
			string country_j = iso3[row_j];
			int region_j = region[row_j];
			int super_region_j = super_region[row_j];



			int this_relatedness ;
			double this_zeta_weight ;
			if (country_i == country_j) {
				this_relatedness = 0;
				this_zeta_weight = this_zeta;
			}
			else if(region_i == region_j) {
				this_relatedness = 1 ;
				this_zeta_weight = this_zeta*(1-this_zeta) ;
			}
			else if(super_region_i == super_region_j) {
				this_relatedness = 2 ;
				this_zeta_weight = pow(1-this_zeta,2) ;
			}
			else {
				this_relatedness = 3;
				this_zeta_weight = 0; // no weight
			}

			double agew,timew,age_time ;


			agew = age_weights.GetWeight(age_i, age_j) ;
			timew = year_weights.GetWeight(year_i, year_j, this_lambda) ;

			age_time = timew*agew ;
			w_mat[j].wt = this_zeta_weight*age_time ;
			w_mat[j].relatedness =(this_relatedness) ;


			row_sums[this_relatedness] += age_time ;
			total_sum += age_time ;
			if (national[row_j] == 0 ) {
				subnational_sum += age_time ;
			}
		}



		for (int j = 0; j < has_data_rows.size() ; ++j) {
			int row_j = has_data_rows[j];

			if (national[row_j] == 0 && fabs(subnational_sum-total_sum) > AGE_EPS) {
				w_mat[j].wt =w_mat[j].wt/subnational_sum * ( (total_sum - subnational_sum) * (1/this_zeta) - (total_sum - subnational_sum) ) ;
			}

			w_mat[j].wt = w_mat[j].wt / row_sums[w_mat[j].relatedness] ;
			if (has_data[row_j]) {
				#pragma omp atomic
				residual_weighted_average_totals[i] += w_mat[j].wt;
				for (int k = 0; k < num_submodels ; ++k) {
					pred_residual[row_num*num_submodels+k] += res_data[row_j*num_submodels+k]*w_mat[j].wt ;
				}
			}
			if (fabs(ln_rate[row_j]) > AGE_EPS) {
				#pragma omp atomic
				simple_weighted_average_totals[i] += w_mat[j].wt;
				simple_spacetime_ln_rate[row_num] += ln_rate[row_j]*w_mat[j].wt ;
			}


		}



	}
	cout << "I made it " << endl;

	double sum_of_simple =std::accumulate(simple_weighted_average_totals.begin(),simple_weighted_average_totals.end(),0.0);
	cout << "Sum simple weighted average totals: " << sum_of_simple;
	double sum_of_residual =std::accumulate(residual_weighted_average_totals.begin(),residual_weighted_average_totals.end(),0.0);
	cout << "Sum residual weighted average totals: " << sum_of_residual;
	#pragma omp parallel for
		for (int i = 0; i < prediction_rows.size() ; ++i) {
			int row_num = prediction_rows[i];
			for (int k = 0; k < num_submodels ; ++k) {
				pred_residual[row_num*num_submodels+k] =pred_residual[row_num*num_submodels+k]/residual_weighted_average_totals[i] ;
			}
			simple_spacetime_ln_rate[row_num] = simple_spacetime_ln_rate[row_num]/simple_weighted_average_totals[i] ;
		}


	//for (int i = 0; i < prediction_rows.size(); ++i ) {
	//	cout << "Residual weighted average" << residual_weighted_average_totals[i] << endl;
	//}

	// write file
	ofstream out;
	out.precision(12) ;

	out.open(outname.c_str());
	cout << "HHi";
	// write the header
	ostringstream temp;
	temp.precision(12);
	temp << "simple_spacetime_ln_rate," ;
	for (int i = 0; i < num_submodels; ++i) {
		temp << "pred_residual_" << i+1  << "," ;
	}
	string line = temp.str() ;
	out << line << "\n" ;

	// write the rows
	for (int i = 0; i < num_obs; ++i) {
		ostringstream temp;
		// matches stata precision
		temp.precision(12);

		// Figure out what the country-age is
		country_age current_pair;
		current_pair.iso3 = iso3[i];
		current_pair.age = age[i];

		// Missing if it should be missing
		if (std::find(prediction_pairs.begin(), prediction_pairs.end(), current_pair)==prediction_pairs.end()) {
			temp << "." << "," ;
		}
		else  {
			temp << simple_spacetime_ln_rate[i] << "," ;
		}

		for (int j = 0; j < num_submodels; ++j) {

			if (std::find(prediction_pairs.begin(), prediction_pairs.end(), current_pair)==prediction_pairs.end()) {
				temp << "." << "," ;
			}
			else if (has_data_total >= 5) {
			temp << pred_residual[i*num_submodels+j] << "," ;
			}
			else {
				temp << "0" << "," ;
			}
		}
		string line = temp.str() ;
		out << line << "\n";
	}

	// close the file
	out.close();
	return 0;
}
