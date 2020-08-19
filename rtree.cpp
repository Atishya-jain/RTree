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
	for(int j = 0; j<nodesize; j++){
		memcpy((data + j*jump_size), (point+j), jump_size);
	}
}

void copy_to_memory(char* data, int* point){
	for (int i = 0; i < nodesize; ++i){
		memcpy((point+i), (data+i*jump_size), sizeof(int));
	}	
}

bool is_contained(int *mbr, int child, vector <int> *points){
	int start;
	if(child == -1){
		start = 2;
	}else{
		start = 2 + 2*dim + 1 + child*(2*dim+1);
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
		start = 2 + 2*dim + 1 + child*(2*dim+1);
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
			return false;
		}
	}
	return true;
}


// // Children inside a node
// struct Child{
// 	int id;
// 	vector< pair<int,int>> MBR;	
// 	Child() : id(-1)
// 	{
// 		for (int i = 0; i < dim; i++){
// 			MBR.push_back(make_pair(INT_MAX, INT_MIN));
// 		}		
// 	}
// };

// // node of the tree
// struct Node{
// 	int id;
// 	int pid;
// 	vector< pair<int,int> > MBR;
// 	// vector< pair<int, vector<pair<int, int>>> > children;
// 	vector<Child> children;
// 	Node() : id(-1), pid(-1)
// 	{
// 		for (int i = 0; i < dim; ++i){
// 			MBR.push_back(make_pair(INT_MAX, INT_MIN));
// 		}
// 		for (int i = 0; i < maxCap; ++i){
// 			// vector<pair<int, int>> temp;
// 			// for (int i = 0; i < dim; i++){
// 			// 	temp.push_back(make_pair(INT_MAX, INT_MIN));
// 			// }
// 			Child child = Child();
// 			// children.push_back(make_pair(-1, temp));
// 			// Child child = Child(dim);
// 			children.push_back(child);
// 		}
// 	}
// };

void print_point(vector<int> *points){
	cout << "Points: ";
	for(int i = 0; i<points->size(); i++){
		cout << points->at(i) << " ";
	}
	cout << endl;
}

void print_mbr(int *mbr, int child){
	// cout << "point: ";
	// for (int i = 0; i < nodesize; ++i)
	// {
	// 	cout << *(mbr+i) << " ";
	// }
	int start;
	if(child == -1){
		start = 2;
	}else{
		start = 2 + 2*dim + 1 + child*(2*dim+1);
	}
	cout << "Min Point: ";
	for(int i = start; i<start+dim; i++){
		cout << *(mbr+i) << " ";
	}
	cout << endl;

	cout << "Max Point: ";
	for(int i = start+dim; i<start+2*dim; i++){
		cout << *(mbr+i) << " ";
	}
	cout << endl;
}

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
	*(point) = INT_MIN;
	*(point+1) = INT_MIN;
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
	cout << "Recursion: " << start << " " << end << endl;
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
			int point[nodesize];
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
				int child_point[nodesize];
				initialise(child_point);
				// Node child_point = Node();
				// Node *child_point = (Node*) data_child;
				copy_to_memory(data_child, child_point);
				// memcpy(&child_point, data_child, sizeof(Node));
				// Setting parent's id and self child MBR
				cout << "Child: " << child_point[0] << endl;
				print_mbr(child_point, -1);
				child_point[1] = point[0];

				// This assigns MBR and id of children into current parent node
				// point.children[counter].MBR = child_point.MBR;
				// point.children[counter].id = child_point.id;

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
				counter++;
				start++;
				// If maxCap children done, break
				if(counter >= maxCap){
					break;
				}
			}

			// if(counter < maxCap){
			// 	int start_child = begin_children+counter*data_per_child;
			// 	for (int i = start_child; i < nodesize; i++){
			// 		point[i] = INT_MIN;
			// 	}
			// 	for(int i = counter; i< maxCap; i++){
			// 		start_child = begin_children+i*data_per_child;
			// 		point[start_child] = -1;
			// 		for(int j = start_child+1; j<=start_child+dim; j++){
			// 			point[i] = INT_MAX;
			// 		}
			// 	}
			// }

			cout << "Page: " << point[0] << endl;
			print_mbr(point, -1);
			copy_to_disk(data, point);
			// memcpy(&data[0], &point, sizeof(Node));
			fh->MarkDirty(point[0]);
			fh->UnpinPage(point[0]);
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

void foo(int ttemp, FileHandler *fhout){
	PageHandler tempph = fhout->PageAt(ttemp);
	char* tempdata = tempph.GetData();
	// Node tempchild = Node();
	// Node *tempchild = (Node*) tempdata;
	// cout << sizeof(Node) << endl;
	// Node *tempchild = new Node();
	int tempchild[nodesize];
	// memcpy(tempchild, &tempdata[0], nodesize*sizeof(int));
	// for (int i = 0; i < nodesize; ++i){
	// 	memcpy(&tempchild[i], &tempdata[i*jump_size], sizeof(int));
	// }
	copy_to_memory(tempdata, tempchild);
	cout << "TEMP: " << tempchild[0] << endl;
	print_mbr(tempchild, -1);
	// free(tempchild);
	fhout->UnpinPage(ttemp);
	fhout->FlushPage(ttemp);
	// cout <<"END\n";
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
	cout << "Starting with BulkLoad\n";

	while(true){
		int point[nodesize];
		initialise(point);

		phout = fhout->NewPage();
		char* dataout = phout.GetData();
		point[0] = phout.GetPageNum();
		point[1] = -1;
		cout << "New point assigned\n";
		int num_children = 0;
		int index = 0;
		while(num_children<maxCap && i<num_points){
			// Check if page left is enough to contain next d dim point
			if (start + jump_size*dim > PAGE_CONTENT_SIZE){
				int page_number= ph.GetPageNum();
				fh.UnpinPage(page_number);
				// Here ideally we shouldn't get any error. Think if error handling needs to be done ----------------------------
				ph=fh.NextPage(page_number);
				data=ph.GetData();
				start = 0;
			}

			// Now from whatever page and start offset, read the point.
			// Here it is guaranteed for the point to be in page
			vector<int> temp;
			cout << "loc ";
			for (int j = start; j<start+jump_size*dim; j+=jump_size){
				int loc;
				memcpy(&loc, &data[j], jump_size);
				temp.push_back(loc);
				cout << loc << " ";
			}
			cout << endl;

			// Assign MBR as a point sized rectangle
			index = 2 + 2*dim + num_children*data_per_child;
			point[index] = -1;
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

		cout << "Page: " << point[0] << endl;
		print_mbr(point, -1);
		// // Here assigning minimum MBR points as INT_MIN doesn't matter bcoz at leaves we don't check MBR of children
		copy_to_disk(dataout, point);

		fhout->MarkDirty(point[0]);
		fhout->UnpinPage(point[0]);
		if(i >= num_points){
			break;
		}
		foo(0, fhout);
	}

	// Unpin the last page that was being used
	fh.UnpinPage(ph.GetPageNum());

	foo(0, fhout);
	// cout << "Starting with Recursion\n";
	recurse_until_root(fhout, start_pg_num, end_pg_num);	
	return ;
}

// void insert(vector <int> *points, FileHandler *fhout, int pg_num){
// 	PageHandler ph = fhout->PageAt(pg_num);
// 	char *data = ph.GetData();
// 	int node[nodesize];
// 	copy_to_memory(data, node);
// 	// If the node is not a leaf
// 	if(!is_leaf(node)){
// 		// Finsing child with min area of expansion on insert
// 		double exp_area;
// 		double min_area = -1;
// 		int min_area_pg = -1;
// 		for (int i = 0; i < maxCap; i++){
// 			exp_area = get_expanded_area(points, node, i);
// 			if (exp_area <= min_area || min_area == -1){
// 				min_area = exp_area;
// 				min_area_pg = get_child_pg(node, i);
// 			}
// 		}
// 		fhout->UnpinPage(pg_num);
// 		// Recurse on that child
// 		insert(points, fhout, min_area_pg);
// 	}else{
// 		// If the node is leaf
// 		// Find number of children already there
// 		int num_children = get_num_children(node);
// 		// If num children less than max, assign, update MBR and return, else split
// 		if (num_children < maxCap){
// 			node[begin_children+num_children*data_per_child] = 
// 		}else{

// 		}
// 	}
// }


bool query(vector <int> *points, FileHandler *fhout, queue<int>*bfs_queue){
	PageHandler ph;
// 	// char* data;
// 	// Node *node;
	int pg_num;

	print_point(points);

	while(!bfs_queue->empty()){
		pg_num = bfs_queue->front();
		bfs_queue->pop();
		cout << "Node ID: " << pg_num << endl;
		ph = fhout->PageAt(pg_num);
		char *data = ph.GetData();
		int node[nodesize];
		copy_to_memory(data, node);
		print_mbr(node, -1);
		if(!is_leaf(node)){
			cout << "Internal\n";
			for(int i = 0; i<maxCap; i++){
				print_mbr(node, i);
				if(is_contained(node, i, points)){
					cout << "Something new in BFS\n";
					bfs_queue->push(node[begin_children+i*(2*dim+1)]);
				}
			}
		}else{
			cout << "Leaf\n";
			for(int i = 0; i<maxCap; i++){
				if(is_contained(node, i, points)){
					cout << "Found\n";
					fhout->UnpinPage(pg_num);
					return true;
				}
			}
		}
		fhout->UnpinPage(pg_num);
		// fm.PrintBuffer();
		cout << "queue size: " << bfs_queue->size() << endl;
	}
	return false;
}

int main(int argc, char *argv[]) {
	string inputfile = argv[1];
	maxCap = stoi(argv[2]);
	dim = stoi(argv[3]);
	string outputfile = argv[4];
	cout << inputfile << " " << maxCap << " " << dim << " " << outputfile << endl;

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
				bulkload(commands[1], stoi(commands[2]), &fhout);
				outfile << "BULKLOAD\n\n\n";
			}else if(commands[0] == "INSERT"){
				vector<int> points;
				for (int i = 0; i < dim; i++){
					points.push_back(stoi(commands[i+1]));
				}
				// insert(&points, &fhout);
				outfile << "INSERT\n\n\n";
			}else{
				vector<int> points;
				for (int i = 0; i < dim; i++){
					points.push_back(stoi(commands[i+1]));
				}
				queue <int> bfs_queue;
				bfs_queue.push(root_page_number);
				bool found = query(&points, &fhout, &bfs_queue);				
				if(found == true){
					outfile << "TRUE\n\n\n";
				}else{
					outfile << "FALSE\n\n\n";
				}
			}
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