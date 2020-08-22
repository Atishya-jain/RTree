#include<iostream>
#include "file_manager.h"
#include "errors.h"
#include "constants.h"
#include<bits/stdc++.h>

using namespace std;

// Globals
// struct Node global_tree;
FileManager fm;
int jump_size = sizeof(int);
int dim, maxCap, nodesize;
int root_page_number = INT_MIN;
int begin_children;
int data_per_child;

void copy_to_disk(char* data, int* point){
	// cout << "Copy to disk starts now\n";
	for(int j = 0; j<nodesize; j++){
		memcpy((data + j*jump_size), (point+j), jump_size);
	}
	// cout << "Copy to disk ends starts now\n";
}

void copy_to_memory(char* data, int* point){
	for (int i = 0; i < nodesize; ++i){
		memcpy((point+i), (data+i*jump_size), sizeof(int));
	}	
}

int get_child_pg(int *node, int child){
	int start;
	if(child == -1){
		start = 2;
	}else{
		start = begin_children + child*data_per_child;
	}
	return *(node+start);	
}

bool is_contained(int *mbr, int child, vector <int> *points){
	int start;
	if(child == -1){
		start = 2;
	}else{
		start = begin_children + 1 + child*data_per_child;
	}

	for(int i = start; i<start+dim; i++){
		if((points->at(i-start) > *(mbr+i+dim)) || (points->at(i-start) < *(mbr+i))){
			return false;
		}
	}
	return true;
}

double get_area(int* node, int child){
	int start;
	if(child == -1){
		start = 2;
	}else{
		start = begin_children + 1 + child*data_per_child;
	}
	double area = 1;
	for(int i = start; i<start+dim; i++){
		area *= (*(node+i+dim)-*(node+i));
	}
	return area;
}

bool is_leaf(int *node){
	for(int i = 0; i<maxCap; i+=1){
		if(get_area(node, i) != 0){
			int pg = get_child_pg(node, i);
			if ( (pg != -1) && (pg != -2) ) {
				return false;
			}

		}
	}
	return true;
}

bool eqDouble(double a, double b){
	if (fabs(a - b) < 0.000000001){
		return true;
	}else{
		return false;
	}
}

pair <double, double> get_expanded_area(vector <int> *points, int *node, int child){
	int start;
	if(child == -1){
		start = 2;
	}else{
		start = begin_children + 1 + child*data_per_child;
	}
	double area = 1;
	double bigger_area = 1;
	for(int i = start; i<start+dim; i++){
		area *= (*(node+i+dim)-*(node+i));
		bigger_area *= (max(points->at(i-start), *(node+i+dim)) - min(points->at(i-start), *(node+i)));
	}
	return make_pair(bigger_area-area, area);
}

pair <double, double> get_mbr_expanded_area(int *mbr, int *node, int child){
	int start;
	if(child == -1){
		start = 2;
	}else{
		start = begin_children + 1 + child*data_per_child;
	}
	double area = 1;
	double bigger_area = 1;
	for(int i = start; i<start+dim; i++){
		area *= (*(node+i+dim)-*(node+i));
		bigger_area *= (max(*(mbr+2+i-start+dim), *(node+i+dim)) - min(*(mbr+2+i-start), *(node+i)));
	}
	return make_pair(bigger_area-area, area);
}

int get_num_children(int *node){
	int num_children = 0;
	int index = begin_children;
	for(int i = begin_children; i<nodesize; i += data_per_child){
		if(*(node+i) == -1){
			return num_children;
		}
		num_children++;
	}
	return num_children;
}

// void print_point(vector<int> *points){
// 	cout << "Points: ";
// 	for(int i = 0; i<points->size(); i++){
// 		cout << points->at(i) << " ";
// 	}
// 	cout << endl;
// }

// void print_mbr(int *mbr, int child){
// 	int start;
// 	if(child == -1){
// 		start = 2;
// 	}else{
// 		start = 2 + 2*dim + 1 + child*(2*dim+1);
// 	}
// 	cout << "Min Point: ";
// 	for(int i = start; i<start+dim; i++){
// 		cout << *(mbr+i) << " ";
// 	}
// 	cout << endl;

// 	cout << "Max Point: ";
// 	for(int i = start+dim; i<start+2*dim; i++){
// 		cout << *(mbr+i) << " ";
// 	}
// 	cout << endl;
// }

void split_line(string line, vector <string> *out){
	string word = "";
	for (int i = 0; i < line.length(); ++i){
		if(line[i] == ' '){
			if(word.length() > 0){
				out->push_back(word);
				word = "";
			}
		}else{
			word += line[i];
		}
	}
	if (word.length() > 0){
		out->push_back(word);		
	}
}

// void getMBRfrompoint(vector<int> point, vector< pair<int,int> > *mbr){
// 	for(int i = 0; i<point.size(); i++){
// 		mbr->at(i) = make_pair(point[i], point[i]);
// 	}
// }

void initialise(int *point){
	*(point) = -1;
	*(point+1) = -1;
	for (int i = 2; i < dim+2; ++i){
		*(point+i) = INT_MAX;
	}
	for (int i = dim+2; i < 2*dim+2; ++i){
		*(point+i) = INT_MIN;
	}
	int start = begin_children;
	for(int i = 0; i<maxCap; i++){
		start = begin_children+i*data_per_child;
		*(point+start) = -1;
		for(int j = start+1; j<=start+dim; j++){
			*(point+j) = INT_MAX;
		}
		for(int j = start+1+dim; j<=start+2*dim; j++){
			*(point+j) = INT_MIN;
		}
	}
}

void recurse_until_root(FileHandler *fh, int start, int end){
	// cout << "Recursion: " << start << " " << end << endl;
	if(end > start){
		PageHandler ph;
		PageHandler ph_child;
		int start_pg_num = INT_MAX;
		int end_pg_num = INT_MIN;
		// char* data;
		// char* data_child;
		int counter = 0;
		while(true){
			ph = fh->NewPage();
			char* data = ph.GetData();
			int *point = new int[nodesize];
			initialise(point);
			// Node point = Node();
			point[0] = ph.GetPageNum();
			// getMBRChild(temp, &point.MBR);
			// Maintain start page as minimum page number 
			if(point[0] <= start_pg_num){
				start_pg_num = point[0];
			}

			// Maintain end page as max page number 
			if(point[0] >= end_pg_num){
				end_pg_num = point[0];
			}

			// Assign all children and set their parent id also
			counter = 0;
			while(start <= end){
				// Fetch this child page and copy into memory
				ph_child = fh->PageAt(start);
				char* data_child = ph_child.GetData();
				int *child_point = new int[nodesize];
				copy_to_memory(data_child, child_point);

				// Setting parent's id and self child MBR
				// cout << "Child: " << child_point[0] << endl;
				// print_mbr(child_point, -1);
				child_point[1] = point[0];

				// This assigns MBR and id of children into current parent node
				int start_child = begin_children+counter*data_per_child;
				point[begin_children+counter*data_per_child] = child_point[0];
				for (int i = start_child+1; i < start_child+data_per_child; i++){
					point[i] = child_point[i - start_child + 1];
				}

				// Update your own MBR
				for(int i = 2; i<dim+2; i++){
					if(point[0] == 100){
						cout << point[i] << " " << child_point[i] << endl;
					}
					point[i] = min(point[i], child_point[i]);
					point[i+dim] = max(point[i+dim], child_point[i+dim]);
				}

				// Copy back child_point into page and mark it dirty
				// memcpy(&data_child[0], &child_point, sizeof(Node));
				copy_to_disk(data_child, child_point);
				fh->MarkDirty(child_point[0]);
				fh->UnpinPage(child_point[0]);
				fh->FlushPage(child_point[0]);
				delete child_point;
				counter++;
				start++;
				// If maxCap children done, break
				if(counter >= maxCap){
					break;
				}
			}

			// cout << "Page: " << point[0] << endl;
			// print_mbr(point, -1);
			copy_to_disk(data, point);
			fh->MarkDirty(point[0]);
			fh->UnpinPage(point[0]);
			fh->FlushPage(point[0]);
			delete point;
			if(start > end){
				break;
			}
		}

		if (start_pg_num == end_pg_num){
			root_page_number = start_pg_num;
		}
		return recurse_until_root(fh, start_pg_num, end_pg_num);
	}
}

void bulkload(string inputfile, int num_points, FileHandler *fhout){
	// Open the input file
	FileHandler fh = fm.OpenFile(inputfile.c_str());
	PageHandler phout;
	// char *dataout;

	// Get first page
	PageHandler ph = fh.FirstPage();
	char *data = ph.GetData ();
	
	// Start your reading
	int start = 0;
	int start_pg_num = INT_MAX;
	int end_pg_num = INT_MIN;
	int i = 0;
	// cout << "Starting with BulkLoad\n";

	while(true){
		int *point = new int[nodesize];
		initialise(point);

		phout = fhout->NewPage();
		char* dataout = phout.GetData();
		point[0] = phout.GetPageNum();
		point[1] = -1;
		// cout << "New point assigned\n";
		int num_children = 0;
		int index = 0;
		while(num_children<maxCap && i<num_points){
			// Check if page left is enough to contain next d dim point
			if (start + jump_size*dim > PAGE_CONTENT_SIZE){
				int page_number= ph.GetPageNum();
				fh.UnpinPage(page_number);
				fh.FlushPage(page_number);
				// Here ideally we shouldn't get any error. Think if error handling needs to be done ----------------------------
				ph=fh.NextPage(page_number);
				data=ph.GetData();
				start = 0;
			}

			// Now from whatever page and start offset, read the point.
			// Here it is guaranteed for the point to be in page
			vector<int> temp;
			// cout << "loc ";
			for (int j = start; j<start+jump_size*dim; j+=jump_size){
				int loc;
				memcpy(&loc, &data[j], jump_size);
				temp.push_back(loc);
				// cout << loc << " ";
			}
			// cout << endl;

			// Assign MBR as a point sized rectangle
			index = begin_children + num_children*data_per_child;
			// -2 ID for special child leaves
			point[index] = -2;
			for(int j = index+1; j<index+1+dim; j++){
				point[j] = temp[j-index-1];
				point[j+dim] = temp[j-index-1];
				point[j-index+1] = min(point[j-index+1], point[j]);
				point[j-index+1+dim] = max(point[j-index+1+dim], point[j+dim]);
			}
			num_children++;
			i++;
			start += jump_size*dim;
		}
		// Maintain start page as minimum page number 
		if(point[0] <= start_pg_num){
			start_pg_num = point[0];
		}

		// Maintain end page as max page number 
		if(point[0] >= end_pg_num){
			end_pg_num = point[0];
		}

		// cout << "Page: " << point[0] << endl;
		// print_mbr(point, -1);
		// Here assigning minimum MBR points as INT_MIN doesn't matter bcoz at leaves we don't check MBR of children
		copy_to_disk(dataout, point);

		fhout->MarkDirty(point[0]);
		fhout->UnpinPage(point[0]);
		fhout->FlushPage(point[0]);
		delete point;
		if(i >= num_points){
			break;
		}
	}

	// Unpin the last page that was being used
	fh.UnpinPage(ph.GetPageNum());
	fh.FlushPage(ph.GetPageNum());

	recurse_until_root(fhout, start_pg_num, end_pg_num);	
	return ;
}

int get_preference(int *old_node, int *new_node, int *node, int child, int child_index, int num_old_children, int num_new_children){
	// If a node requires all remaining children to satisfy min capacity, assign
	int min_reqd = ceil(maxCap/2);
	// maxCap-child-2 number of children left to allocate including this
	if((min_reqd - num_old_children) >= (maxCap-child-1)){
		return 0;
	}
	if((min_reqd - num_new_children) >= (maxCap-child-1)){
		return 1;
	}

	pair<double, double> old_area = get_mbr_expanded_area(old_node, node, child_index);
	pair<double, double> new_area = get_mbr_expanded_area(new_node, node, child_index);

	if(old_area.first < new_area.first){
		return 0;
	}else if (eqDouble(old_area.first, new_area.first)){
		if(old_area.second <new_area.second){
			return 0;
		}else if(eqDouble(old_area.second, new_area.second)){
			if(num_old_children < num_new_children){
				return 0;
			}else{
				return 1;
			}
		}else{
			return 1;
		}
	}else{
		return 1;
	}
}

int get_preference_points(int *old_node, int *new_node, vector<int> *points, int num_old_children, int num_new_children){
	// If a node requires all remaining children to satisfy min capacity, assign
	int min_reqd = ceil(maxCap/2);
	// maxCap-child-2 number of children left to allocate including this
	if((min_reqd - num_old_children) >= 1){
		return 0;
	}
	if((min_reqd - num_new_children) >= 1){
		return 1;
	}

	pair<double, double> old_area = get_expanded_area(points, old_node, -1);
	pair<double, double> new_area = get_expanded_area(points, new_node, -1);

	if(old_area.first < new_area.first){
		return 0;
	}else if (eqDouble(old_area.first, new_area.first)){
		if(old_area.second <new_area.second){
			return 0;
		}else if(eqDouble(old_area.second, new_area.second)){
			if(num_old_children < num_new_children){
				return 0;
			}else{
				return 1;
			}
		}else{
			return 1;
		}
	}else{
		return 1;
	}
}

double get_mbr_area(int *node, int i, int j){
	double bigarea = 1;
	double area1 = 1;
	double area2 = 1;
	int index1 = begin_children+i*data_per_child;
	int index2 = begin_children+j*data_per_child;

	for(int k = 0; k<dim; k++){
		area1 *= (*(node+index1+dim+k) - *(node+index1+k));
		area2 *= (*(node+index2+dim+k) - *(node+index2+k));
		bigarea *= (max(*(node+index1+dim+k), *(node+index2+dim+k)) - min(*(node+index1+k), *(node+index2+k)));
	}
	return bigarea - area1 - area2;
}

double get_mbr_area3(int *node, int i, int *mbr){
	double bigarea = 1;
	double area1 = 1;
	double area2 = 1;
	int index1 = begin_children+i*data_per_child;
	int index2 = 2;

	for(int k = 0; k<dim; k++){
		area1 *= (*(node+index1+dim+k) - *(node+index1+k));
		area2 *= (*(node+index2+dim+k) - *(node+index2+k));
		bigarea *= (max(*(node+index1+dim+k), *(node+index2+dim+k)) - min(*(node+index1+k), *(node+index2+k)));
	}
	return bigarea - area1 - area2;
}

double get_mbr_area2(int *node, int i, vector<int> *points){
	double bigarea = 1;
	int index1 = begin_children+i*data_per_child;
	int index2 = 0;

	for(int k = 0; k<dim; k++){
		bigarea *= (max(*(node+index1+dim+k), points->at(k)) - min(*(node+index1+k), points->at(k)));
	}
	return bigarea;
}


pair<int, int> get_seeds2(int *node, int *newnode){
	double max_dist = INT_MIN;
	double running_dist = INT_MIN;
	pair<int, int> answer = make_pair(-1,-1);
	for (int i = 0; i < maxCap; i++){
		for(int j = i+1; j <= maxCap; j++){
			if(j != maxCap){
				running_dist = get_mbr_area(node, i, j);
				if (running_dist > max_dist){
					max_dist = running_dist;
					answer.first = i;
					answer.second = j;
				}
			}else{
				running_dist = get_mbr_area3(node, i, newnode);
				if (running_dist > max_dist){
					max_dist = running_dist;
					answer.first = i;
					answer.second = j;
				}
			}
		}
	}
	return answer;
}

pair<int, int> get_seeds(int *node, vector<int> *points){
	double max_dist = INT_MIN;
	double running_dist = INT_MIN;
	pair<int, int> answer = make_pair(-1,-1);
	for (int i = 0; i < maxCap; i++){
		for(int j = i+1; j <= maxCap; j++){
			if(j != maxCap){
				running_dist = get_mbr_area(node, i, j);
				if (running_dist > max_dist){
					max_dist = running_dist;
					answer.first = i;
					answer.second = j;
				}
			}else{
				running_dist = get_mbr_area2(node, i, points);
				if (running_dist > max_dist){
					max_dist = running_dist;
					answer.first = i;
					answer.second = j;
				}
			}
		}
	}
	return answer;
}

void AdjustParent(int old_pg_num, int new_pg_num, FileHandler *fhout){
	PageHandler old_ph = fhout->PageAt(old_pg_num);
	char *old_data = old_ph.GetData();
	int *old_node = new int[nodesize];
	copy_to_memory(old_data, old_node);

	PageHandler new_ph = fhout->PageAt(new_pg_num);
	char *new_data = new_ph.GetData();
	int *new_node = new int[nodesize];
	copy_to_memory(new_data, new_node);

	int par_pg_num = new_node[1];
	// cout << "Is it Root: " << old_pg_num << " " << new_pg_num << " " << par_pg_num << endl;
	// print_mbr(old_node, -1);
	// print_mbr(new_node, -1);

	if (par_pg_num == -1){
		PageHandler root = fhout->NewPage();
		char *root_data = root.GetData();
		int root_pg_num = root.GetPageNum();
		root_page_number = root_pg_num;
		int *root_node = new int[nodesize];
		initialise(root_node);
		// copy_to_memory(root_data, root_node);

		old_node[1] = root_pg_num;
		new_node[1] = root_pg_num;

		root_node[0] = root_page_number;
		root_node[1] = -1;
		int index = begin_children;
		root_node[index] = old_node[0];
		for(int j = index+1; j<index+1+dim; j++){
			root_node[j] = old_node[j-index+1];
			root_node[j+dim] = old_node[j-index+1+dim];
			root_node[j-index+1] = min(root_node[j-index+1], old_node[j-index+1]);
			root_node[j-index+1+dim] = max(root_node[j-index+1+dim], old_node[j-index+1+dim]);
		}

		index += data_per_child; 
		root_node[index] = new_node[0];
		for(int j = index+1; j<index+1+dim; j++){
			root_node[j] = new_node[j-index+1];
			root_node[j+dim] = new_node[j-index+1+dim];
			root_node[j-index+1] = min(root_node[j-index+1], new_node[j-index+1]);
			root_node[j-index+1+dim] = max(root_node[j-index+1+dim], new_node[j-index+1+dim]);
		}

		// cout << "Root MBR: " << root_node[begin_children] << " " << root_node[begin_children+data_per_child] << endl;
		// print_mbr(root_node, -1);


		copy_to_disk(new_data, new_node);
		copy_to_disk(old_data, old_node);
		copy_to_disk(root_data, root_node);
		fhout->MarkDirty(new_pg_num);
		fhout->MarkDirty(old_pg_num);
		fhout->MarkDirty(root_pg_num);
		fhout->UnpinPage(old_pg_num);
		fhout->UnpinPage(new_pg_num);
		fhout->UnpinPage(root_pg_num);
		fhout->FlushPage(old_pg_num);
		fhout->FlushPage(new_pg_num);
		fhout->FlushPage(root_pg_num);
		delete new_node;
		delete old_node;
		delete root_node;
	}else{
		PageHandler par_ph = fhout->PageAt(par_pg_num);
		char *par_data = par_ph.GetData();		
		int *par_node = new int[nodesize];
		copy_to_memory(par_data, par_node);
		int num_children = get_num_children(par_node);

		// Correct child's credentials here in parent
		for(int i = 0; i<num_children; i++){
			if(get_child_pg(par_node, i) == old_pg_num){
				int index = begin_children + i*data_per_child;
				for(int j = index+1; j<index+1+dim; j++){
					par_node[j] = old_node[j-index+1];
					par_node[j+dim] = old_node[j-index+1+dim];
					par_node[j-index+1] = min(par_node[j-index+1], old_node[j-index+1]);
					par_node[j-index+1+dim] = max(par_node[j-index+1+dim], old_node[j-index+1]);
				}
				break;
			}
		}
		// cout << "Page Adjust: " << num_children << endl;

		if (num_children < maxCap){
			int index = begin_children + num_children*data_per_child;
			par_node[index] = new_node[0];
			for(int j = index+1; j<index+1+dim; j++){
				par_node[j] = new_node[j-index+1];
				par_node[j+dim] = new_node[j-index+1+dim];
				par_node[j-index+1] = min(par_node[j-index+1], new_node[j-index+1]);
				par_node[j-index+1+dim] = max(par_node[j-index+1+dim], new_node[j-index+1+dim]);
			}

			// cout << "MBR Adjust: \n";
			// print_mbr(par_node, -1);

			copy_to_disk(new_data, new_node);
			copy_to_disk(par_data, par_node);
			fhout->MarkDirty(new_pg_num);
			fhout->MarkDirty(par_pg_num);
			fhout->UnpinPage(new_pg_num);
			fhout->UnpinPage(par_pg_num);
			fhout->UnpinPage(old_pg_num);
			fhout->FlushPage(new_pg_num);
			fhout->FlushPage(par_pg_num);
			fhout->FlushPage(old_pg_num);
			delete new_node;
			delete old_node;
			delete par_node;
		}else{
			// // Here we would need to read old node as well to redistribute
			// PageHandler old_ph = fhout->PageAt(old_pg_num);
			// char *old_data = old_ph.GetData();
			// int old_node[nodesize];
			// copy_to_memory(old_data, old_node);

			// Create new parent page and two new nodes to be made
			PageHandler new_par_ph = fhout->NewPage();
			char *new_par_data = new_par_ph.GetData();
			int new_par_pg_num = new_par_ph.GetPageNum();

			int *new_par_node = new int[nodesize];
			initialise(new_par_node);

			int *old_par_node = new int[nodesize];
			initialise(old_par_node);

			// Choose two farthest seeds
			pair<int, int> seeds = get_seeds2(par_node, new_node);

			// Form MBR and Children of these two split nodes
			int start_old = 1;
			int start_new = 1;

			// Assigning children and updating MBR simultaneously
			old_par_node[0] = par_node[0];
			old_par_node[1] = par_node[1];
			int index = begin_children + seeds.first*data_per_child;
			old_par_node[begin_children] = par_node[index];
			for(int j = index+1; j<index+1+dim; j++){
				old_par_node[begin_children+j-index] = par_node[j];
				old_par_node[begin_children+j-index+dim] = par_node[j+dim];
				old_par_node[j-index+1] = min(old_par_node[j-index+1], par_node[j]);
				old_par_node[j-index+1+dim] = max(old_par_node[j-index+1+dim], par_node[j+dim]);
			}

			new_par_node[0] = new_par_pg_num;
			new_par_node[1] = par_node[1];
			if(seeds.second < maxCap){
				index = begin_children + seeds.second*data_per_child;
				new_par_node[begin_children] = par_node[index];
				for(int j = index+1; j<index+1+dim; j++){
					new_par_node[begin_children+j-index] = par_node[j];
					new_par_node[begin_children+j-index+dim] = par_node[j+dim];
					new_par_node[j-index+1] = min(new_par_node[j-index+1], par_node[j]);
					new_par_node[j-index+1+dim] = max(new_par_node[j-index+1+dim], par_node[j+dim]);
				}
				// cout << "Page Requested of Child: " << par_node[index] << endl;
				PageHandler temph = fhout->PageAt(par_node[index]);
				char* tempdata = temph.GetData();
				int *tempnode = new int[nodesize];
				copy_to_memory(tempdata, tempnode);
				tempnode[1] = new_par_node[0];
				copy_to_disk(tempdata, tempnode);
				fhout->MarkDirty(par_node[index]);
				fhout->UnpinPage(par_node[index]);				
				fhout->FlushPage(par_node[index]);
				delete tempnode;				
			}else{
				// TODO Allocate new_node
				new_par_node[begin_children] = new_node[0];
				for(int j = 1; j<1+dim; j++){
					new_par_node[begin_children+j] = new_node[j+1];
					new_par_node[begin_children+j+dim] = new_node[j+1+dim];
					new_par_node[j+1] = min(new_par_node[j+1], new_node[j+1]);
					new_par_node[j+1+dim] = max(new_par_node[j+1+dim], new_node[j+1+dim]);
				}
				// Update the page number of new_node
				new_node[1] = new_par_node[0];
			}

			int num = -1;
			int index2;
			int num_child = 0;
			for (int i = 0; i < maxCap; i++){
				if(i != seeds.first && i != seeds.second){
					// This also checks if one of the child requires all children to maintain min m children
					num = get_preference(old_par_node, new_par_node, par_node, num_child, i, start_old, start_new);
					// cout << "Considering point: " << get_child_pg(par_node, i) << endl;
					// print_mbr(par_node, i);
					if(num == 0){
						index = begin_children + i*data_per_child;
						index2 = begin_children + start_old*data_per_child;
						old_par_node[index2] = par_node[index];
						for(int j = index+1; j<index+1+dim; j++){
							old_par_node[index2+j-index] = par_node[j];
							old_par_node[index2+j-index+dim] = par_node[j+dim];
							old_par_node[j-index+1] = min(old_par_node[j-index+1], par_node[j]);
							old_par_node[j-index+1+dim] = max(old_par_node[j-index+1+dim], par_node[j+dim]);
						}
						// cout << "Old new MBR: " << endl;
						// print_mbr(old_par_node, -1);
						start_old += 1;
					}else{
						index = begin_children + i*data_per_child;
						index2 = begin_children + start_new*data_per_child;
						new_par_node[index2] = par_node[index];
						for(int j = index+1; j<index+1+dim; j++){
							new_par_node[index2+j-index] = par_node[j];
							new_par_node[index2+j-index+dim] = par_node[j+dim];
							new_par_node[j-index+1] = min(new_par_node[j-index+1], par_node[j]);
							new_par_node[j-index+1+dim] = max(new_par_node[j-index+1+dim], par_node[j+dim]);
						}
						// cout << "New new MBR: " << endl;
						// print_mbr(new_par_node, -1);

						start_new += 1;
						// cout << "Page Requested of Childs: " << par_node[index] << endl;
						PageHandler temph = fhout->PageAt(par_node[index]);
						char* tempdata = temph.GetData();
						int *tempnode = new int[nodesize];
						copy_to_memory(tempdata, tempnode);
						tempnode[1] = new_par_node[0];
						copy_to_disk(tempdata, tempnode);
						fhout->MarkDirty(par_node[index]);
						fhout->UnpinPage(par_node[index]);
						fhout->FlushPage(par_node[index]);
						delete tempnode;
					}
					num_child++;
				}
			}

			// cout << "Kids split Till now: " << start_old << " " << start_new << endl;
			if(seeds.second != maxCap){
				num = get_preference(old_par_node, new_par_node, new_node, num_child, -1, start_old, start_new);
				if(num == 0){
					index2 = begin_children + start_old*data_per_child;
					old_par_node[index2] = new_node[0];
					for(int j = 1; j<1+dim; j++){
						old_par_node[index2+j] = new_node[j+1];
						old_par_node[index2+j+dim] = new_node[j+1+dim];
						old_par_node[j+1] = min(old_par_node[j+1], new_node[j+1]);
						old_par_node[j+1+dim] = max(old_par_node[j+1+dim], new_node[j+1+dim]);
					}
					start_old += 1;
					// Update the page number of new_node
					new_node[1] = old_par_node[0];
				}else{
					index2 = begin_children + start_new*data_per_child;
					new_par_node[index2] = new_node[0];
					for(int j = 1; j<1+dim; j++){
						new_par_node[index2+j] = new_node[j+1];
						new_par_node[index2+j+dim] = new_node[j+1+dim];
						new_par_node[j+1] = min(new_par_node[j+1], new_node[j+1]);
						new_par_node[j+1+dim] = max(new_par_node[j+1+dim], new_node[j+1+dim]);
					}
					start_new += 1;
					new_node[1] = new_par_node[0];
				}
			}

			// cout << "Kids split: " << start_old << " " << start_new << endl;
			// cout << "Old Adjust MBR: " << old_par_node[0] << endl;
			// print_mbr(old_par_node, -1);
			// cout << "New Adjust MBR: " << new_par_node[0] << endl;
			// print_mbr(new_par_node, -1);
			
			// print_mbr(new_par_node);


			copy_to_disk(new_data, new_node);
			copy_to_disk(par_data, old_par_node);
			copy_to_disk(new_par_data, new_par_node);
			fhout->MarkDirty(new_pg_num);
			fhout->MarkDirty(par_pg_num);
			fhout->MarkDirty(new_par_pg_num);
			fhout->UnpinPage(new_pg_num);
			fhout->UnpinPage(par_pg_num);
			fhout->UnpinPage(new_par_pg_num);
			fhout->UnpinPage(old_pg_num);

			fhout->FlushPage(new_pg_num);
			fhout->FlushPage(par_pg_num);
			fhout->FlushPage(new_par_pg_num);
			fhout->FlushPage(old_pg_num);

			delete new_node;
			delete old_node;
			delete par_node;
			delete new_par_node;
			delete old_par_node;
			// Load Parent Page
			AdjustParent(par_pg_num, new_par_pg_num, fhout);
		}
	}
}

void insert(vector <int> *points, FileHandler *fhout, int pg_num){
	// cout << "Page Number: " << pg_num << endl; 
	PageHandler ph = fhout->PageAt(pg_num);
	char *data = ph.GetData();
	int *node = new int[nodesize];
	copy_to_memory(data, node);
	// cout << "My Parent Page Number: " << node[1] << endl; 
	// print_point(points);
	// If the node is not a leaf
	if(!is_leaf(node)){
		// cout << "Non Leaf: " << pg_num << endl;
		// Finsing child with min area of expansion on insert
		pair<double, double> exp_area;
		double min_area = -1;
		double min_area_total;
		int min_area_pg = -1;
		int min_area_index = -1;
		int num_children = get_num_children(node);
		// cout << "num_children: " << num_children << endl;
		for (int i = 0; i < num_children; i++){
			exp_area = get_expanded_area(points, node, i);
			if ((exp_area.first < min_area) || (eqDouble(exp_area.first ,min_area) && exp_area.second < min_area_total) || min_area == -1){
				min_area = exp_area.first;
				min_area_total = exp_area.second;
				min_area_pg = get_child_pg(node, i);
				min_area_index = i;
				// cout << "Exp_Area: " << exp_area.first << " " << exp_area.second << " " << min_area_pg << endl;
			}
		}
		// Update the MBR of Rtree during downward pass as this would become this only in case split doesn't propagate till this stage
		// cout << "Selected Page: " << min_area_pg << endl;
		int index = begin_children + min_area_index*data_per_child;
		for(int j = index+1; j<index+dim+1; j++){
			node[j] = min(node[j], points->at(j-index-1));
			node[j+dim] = max(node[j+dim], points->at(j-index-1));
			node[j-index+1] = min(node[j-index+1], node[j]);
			node[j-index+1+dim] = max(node[j-index+1+dim], node[j+dim]);
		}
		// cout << "Going Down: \n";
		// print_mbr(node, -1);
		// print_mbr(node, min_area_index);

		copy_to_disk(data, node);
		fhout->MarkDirty(pg_num);
		fhout->UnpinPage(pg_num);
		fhout->FlushPage(pg_num);
		delete node;
		// Recurse on that child
		insert(points, fhout, min_area_pg);
	}else{
		// If the node is leaf
		// Find number of children already there
		int num_children = get_num_children(node);
		// cout << "Leaf: " << num_children << endl;
		// If num children less than max, assign, update MBR and return, else split
		if (num_children < maxCap){
			int index = begin_children+num_children*data_per_child;
			node[index] = -2;
			// Assigning children and updating MBR simultaneously
			for(int j = index+1; j<index+1+dim; j++){
				node[j] = points->at(j-index-1);
				node[j+dim] = points->at(j-index-1);
				node[j-index+1] = min(node[j-index+1], node[j]);
				node[j-index+1+dim] = max(node[j-index+1+dim], node[j+dim]);
			}

			// print_mbr(node, -1);
			copy_to_disk(data, node);
			fhout->MarkDirty(pg_num);
			fhout->UnpinPage(pg_num);
			fhout->FlushPage(pg_num);
			delete node;
			// fm.PrintBuffer();
		}else{
			// Create a new page
			PageHandler ph_new = fhout->NewPage();
			char *new_data = ph_new.GetData();
			int new_pg_num = ph_new.GetPageNum();

			int *new_node = new int[nodesize];
			initialise(new_node);

			int *old_node = new int[nodesize];
			initialise(old_node);

			// Choose two farthest seeds
			pair<int, int> seeds = get_seeds(node, points);

			// Form MBR and Children of these two split nodes
			int start_old = 1;
			int start_new = 1;

			// Assigning children and updating MBR simultaneously
			old_node[0] = node[0];
			old_node[1] = node[1];
			int index = begin_children + seeds.first*data_per_child;
			old_node[begin_children] = -2;
			for(int j = index+1; j<index+1+dim; j++){
				old_node[begin_children+j-index] = node[j];
				old_node[begin_children+j-index+dim] = node[j+dim];
				old_node[j-index+1] = min(old_node[j-index+1], node[j]);
				old_node[j-index+1+dim] = max(old_node[j-index+1+dim], node[j+dim]);
			}

			new_node[0] = new_pg_num;
			new_node[1] = node[1];
			if(seeds.second < maxCap){
				index = begin_children + seeds.second*data_per_child;
				new_node[begin_children] = -2;
				for(int j = index+1; j<index+1+dim; j++){
					new_node[begin_children+j-index] = node[j];
					new_node[begin_children+j-index+dim] = node[j+dim];
					new_node[j-index+1] = min(new_node[j-index+1], node[j]);
					new_node[j-index+1+dim] = max(new_node[j-index+1+dim], node[j+dim]);
				}
			}else{
				new_node[begin_children] = -2;
				for(int j = 1; j<=dim; j++){
					new_node[begin_children+j] = points->at(j-1);
					new_node[begin_children+j+dim] = points->at(j-1);
					new_node[j+1] = min(new_node[j+1], points->at(j-1));
					new_node[j+1+dim] = max(new_node[j+1+dim], points->at(j-1));
				}
			}

			int num = -1;
			int index2;
			int num_child = 0;
			for (int i = 0; i < maxCap; i++){
				if(i != seeds.first && i != seeds.second){
					// This also checks if one of the child requires all children to maintain min m children
					num = get_preference(old_node, new_node, node, num_child, i, start_old, start_new);
					// cout << "Considering point: " << get_child_pg(node, i) << endl;
					// print_mbr(node, i);

					if(num == 0){
						index = begin_children + i*data_per_child;
						index2 = begin_children + start_old*data_per_child;
						old_node[index2] = -2;
						for(int j = index+1; j<index+1+dim; j++){
							old_node[index2+j-index] = node[j];
							old_node[index2+j-index+dim] = node[j+dim];
							old_node[j-index+1] = min(old_node[j-index+1], node[j]);
							old_node[j-index+1+dim] = max(old_node[j-index+1+dim], node[j+dim]);
						}
						// cout << "Old new MBR: " << endl;
						// print_mbr(old_node, -1);
						start_old += 1;
					}else{
						index = begin_children + i*data_per_child;
						index2 = begin_children + start_new*data_per_child;
						new_node[index2] = -2;
						for(int j = index+1; j<index+1+dim; j++){
							new_node[index2+j-index] = node[j];
							new_node[index2+j-index+dim] = node[j+dim];
							new_node[j-index+1] = min(new_node[j-index+1], node[j]);
							new_node[j-index+1+dim] = max(new_node[j-index+1+dim], node[j+dim]);
						}
						// cout << "New new MBR: " << endl;
						// print_mbr(new_node, -1);
						start_new += 1;
					}
					num_child++;
				}
			}

			// cout << "Kids split Till now: " << start_old << " " << start_new << endl;
			if(seeds.second != maxCap){
				num = get_preference_points(old_node, new_node, points, start_old, start_new);
				// TODO Points
				if(num == 0){
					index2 = begin_children + start_old*data_per_child;
					old_node[index2] = -2;
					for(int j = 1; j<1+dim; j++){
						old_node[index2+j] = points->at(j-1);
						old_node[index2+j+dim] = points->at(j-1);
						old_node[j+1] = min(old_node[j+1], points->at(j-1));
						old_node[j+1+dim] = max(old_node[j+1+dim], points->at(j-1));
					}
					start_old += 1;
				}else{
					index2 = begin_children + start_new*data_per_child;
					new_node[index2] = -2;
					for(int j = 1; j<1+dim; j++){
						new_node[index2+j] = points->at(j-1);
						new_node[index2+j+dim] = points->at(j-1);
						new_node[j+1] = min(new_node[j+1], points->at(j-1));
						new_node[j+1+dim] = max(new_node[j+1+dim], points->at(j-1));
					}
					start_new += 1;
				}				
			}

			// cout << "Kids split: " << start_old << " " << start_new << endl;
			// cout << "Old MBR: " << old_node[0] << endl;
			// print_mbr(old_node, -1);
			// cout << "New MBR: " << new_node[0] << endl;
			// print_mbr(new_node, -1);

			copy_to_disk(data, old_node);
			copy_to_disk(new_data, new_node);
			fhout->MarkDirty(pg_num);
			fhout->MarkDirty(new_pg_num);
			fhout->UnpinPage(pg_num);
			fhout->UnpinPage(new_pg_num);

			fhout->FlushPage(pg_num);
			fhout->FlushPage(new_pg_num);
			delete node;
			delete old_node;
			delete new_node;
			// Load Parent Page
			AdjustParent(pg_num, new_pg_num, fhout);
		}
	}
}


bool query(vector <int> *points, FileHandler *fhout, queue<int>*bfs_queue){
	PageHandler ph;
// 	// char* data;
// 	// Node *node;
	int pg_num;

	// print_point(points);

	while(!bfs_queue->empty()){
		pg_num = bfs_queue->front();
		bfs_queue->pop();
		// cout << "Node ID: " << pg_num << endl;
		ph = fhout->PageAt(pg_num);
		char *data = ph.GetData();
		int *node = new int[nodesize];
		copy_to_memory(data, node);
		// print_mbr(node, -1);
		if(!is_leaf(node)){
			// cout << "Internal\n";
			int num_children = get_num_children(node);
			for(int i = 0; i<num_children; i++){
				// cout << "Considering child: " << get_child_pg(node, i) << endl;
				// print_mbr(node, i);
				if(is_contained(node, i, points)){
					// cout << "Something new in BFS\n";
					bfs_queue->push(node[begin_children+i*(2*dim+1)]);
				}
			}
		}else{
			// cout << "Leaf\n";
			for(int i = 0; i<maxCap; i++){
				// cout << "Leaf Kid: \n";
				// print_mbr(node, i);
				if(is_contained(node, i, points)){
					// cout << "Found\n";
					fhout->UnpinPage(pg_num);
					fhout->FlushPage(pg_num);
					return true;
				}
			}
		}
		fhout->UnpinPage(pg_num);
		fhout->FlushPage(pg_num);
		delete node;
		// fm.PrintBuffer();
		// cout << "queue size: " << bfs_queue->size() << endl;
	}
	return false;
}

int main(int argc, char *argv[]) {
	string inputfile = argv[1];
	maxCap = stoi(argv[2]);
	dim = stoi(argv[3]);
	string outputfile = argv[4];
	// cout << inputfile << " " << maxCap << " " << dim << " " << outputfile << endl;

	nodesize = 2 + 2*dim + maxCap*(1+2*dim);
	begin_children = 2*dim+2;
	data_per_child = 1+2*dim;

	// global_tree = Node(dim, maxCap);

	ifstream infile(inputfile);
	ofstream outfile (outputfile);
	FileHandler fhout;
	try{
		fhout = fm.CreateFile("bulkload.txt");
	}catch(const InvalidFileException& e){
		fm.DestroyFile("bulkload.txt");
		fhout = fm.CreateFile("bulkload.txt");
	}

	if (infile.is_open() && outfile.is_open()){
		string line;
		while(getline(infile,line)){
			vector <string>commands;
			split_line(line, &commands);
			// cout << commands[0] << " " << commands.size() << '\n';
			if (commands[0] == "BULKLOAD"){
				// cout << "------------------------- BULKLOAD --------------------------------\n";
				bulkload(commands[1], stoi(commands[2]), &fhout);
				outfile << "BULKLOAD\n\n\n";
			}else if(commands[0] == "INSERT"){
				vector<int> points;
				for (int i = 0; i < dim; i++){
					points.push_back(stoi(commands[i+1]));
				}
				// cout << "------------------------- INSERT --------------------------------\n";
				// print_point(&points);
				insert(&points, &fhout, root_page_number);
				outfile << "INSERT\n\n\n";
			}else{
				vector<int> points;
				for (int i = 0; i < dim; i++){
					points.push_back(stoi(commands[i+1]));
				}
				queue <int> bfs_queue;
				bfs_queue.push(root_page_number);
				// cout << "------------------------- QUERY --------------------------------\n";
				bool found = query(&points, &fhout, &bfs_queue);				
				if(found == true){
					outfile << "TRUE\n\n\n";
				}else{
					outfile << "FALSE\n\n\n";
				}
			}
			// fm.PrintBuffer();
		}
		infile.close();
		outfile.close();
	}
	fm.DestroyFile("bulkload.txt");
	// FileManager fm;
	// FileHandler fh = fm.OpenFile("Testcases/TC_6/sortedData_2_10_100.txt");
	// cout << "first page num " << fh.FirstPage().GetPageNum() << endl;
	// fm.CloseFile(fh);

	// fm.PrintBuffer();
	// fm.ClearBuffer();
	// fm.PrintBuffer();

	// FileHandler fh2=fm.OpenFile("Testcases/TC_1/sortedData_2_10_100.txt" );
	// PageHandler ph2;
	// ph2=fh2.NewPage(); //Working fine
	// ph2=fh2.FirstPage(); // gives error: BufferManagerException : Read request failed
}