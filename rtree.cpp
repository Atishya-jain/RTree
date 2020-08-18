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

// void print_point(vector<int> *points){
// 	cout << "Points: ";
// 	for(int i = 0; i<points->size(); i++){
// 		cout << points->at(i) << " ";
// 	}
// 	cout << endl;
// }

void print_mbr(int *mbr){
	// cout << "point: ";
	// for (int i = 0; i < nodesize; ++i)
	// {
	// 	cout << *(mbr+i) << " ";
	// }
	cout << "Max Point: ";
	for(int i = dim+1; i<=2*dim; i++){
		cout << *(mbr+i) << " ";
	}
	cout << endl;

	cout << "Min Point: ";
	for(int i = 1; i<=dim; i++){
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

void getMBRfrompoint(vector<int> point, vector< pair<int,int> > *mbr){
	for(int i = 0; i<point.size(); i++){
		mbr->at(i) = make_pair(point[i], point[i]);
	}
}

// void recurse_until_root(FileHandler *fh, int start, int end){
// 	cout << "Recursion: " << start << " " << end << endl;
// 	if(end > start){
// 		PageHandler ph;
// 		PageHandler ph_child;
// 		int start_pg_num = INT_MAX;
// 		int end_pg_num = INT_MIN;
// 		// char* data;
// 		// char* data_child;
// 		int counter = 0;
// 		while(true){
// 			ph = fh->NewPage();
// 			char* data = ph.GetData();
// 			Node point = Node();
// 			point.id = ph.GetPageNum();
// 			// getMBRChild(temp, &point.MBR);
// 			// Maintain start page as minimum page number 
// 			if(point.id <= start_pg_num){
// 				start_pg_num = point.id;
// 			}

// 			// Maintain end page as max page number 
// 			if(point.id >= end_pg_num){
// 				end_pg_num = point.id;
// 			}

// 			// Assign all children and set their parent id also
// 			counter = 0;
// 			while(start <= end){
// 				// Fetch this child page and copy into memory
// 				ph_child = fh->PageAt(start);
// 				char* data_child = ph_child.GetData();
// 				Node child_point = Node();
// 				// Node *child_point = (Node*) data_child;
// 				memcpy(&child_point, data_child, sizeof(Node));
// 				// Setting parent's id and self child MBR
// 				cout << "Child: " << child_point.id << endl;
// 				print_mbr(&child_point.MBR);
// 				child_point.pid = point.id;
// 				point.children[counter].MBR = child_point.MBR;
// 				point.children[counter].id = child_point.id;
// 				// Update your own MBR
// 				for(int i = 0; i<point.MBR.size(); i++){
// 					point.MBR[i].first = min(point.MBR[i].first, child_point.MBR[i].first);
// 					point.MBR[i].second = max(point.MBR[i].second, child_point.MBR[i].second);
// 				}
// 				// Copy back child_point into page and mark it dirty
// 				memcpy(&data_child[0], &child_point, sizeof(Node));
// 				fh->MarkDirty(child_point.id);
// 				fh->UnpinPage(child_point.id);
// 				counter++;
// 				start++;
// 				// If maxCap children done, break
// 				if(counter >= maxCap){
// 					break;
// 				}
// 			}

// 			cout << "Page: " << point.id << endl;
// 			print_mbr(&point.MBR);

// 			memcpy(&data[0], &point, sizeof(Node));
// 			fh->MarkDirty(point.id);
// 			fh->UnpinPage(point.id);
// 			if(start > end){
// 				break;
// 			}
// 		}

// 		if (start_pg_num == end_pg_num){
// 			root_page_number = start_pg_num;
// 		}
// 		return recurse_until_root(fh, start_pg_num, end_pg_num);
// 	}
// }

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
	print_mbr(tempchild);
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
	int id = 0;
	int start = 0;
	int start_pg_num = INT_MAX;
	int end_pg_num = INT_MIN;

	cout << "Starting with BulkLoad\n";
	for(int i = 0; i<num_points; i+=1){
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
			id += 1;
		}
		cout << endl;
		// Node point = Node();
		int point[nodesize];
		point[0] = id;
		cout << "New point assigned\n";
		// Assign MBR as a point sized rectangle
		for(int j = 1; j<=temp.size(); j++){
			point[j] = temp[j-1];
			point[j+dim] = temp[j-1];
		}

		// getMBRfrompoint(temp, &point.MBR);
		// leaves.push_back(point);
		phout = fhout->NewPage();
		char* dataout = phout.GetData();
		point[0] = phout.GetPageNum();
		// Maintain start page as minimum page number 
		if(point[0] <= start_pg_num){
			start_pg_num = point[0];
		}

		// Maintain end page as max page number 
		if(point[0] >= end_pg_num){
			end_pg_num = point[0];
		}
		cout << "Page: " << point[0] << endl;
		print_mbr(point);

		copy_to_disk(dataout, point);
		// for(int j = 0; j<nodesize; j++){
		// 	memcpy(&dataout[j*jump_size], &point[j], sizeof(int));
		// }
		// memcpy(&dataout[0], &point[0], sizeof(int));
		// int begin = jump_size;
		// for(int j = 0; j<dim; j++){
			// memcpy(&dataout[begin], &point.MBR[j].first, sizeof(int));
			// memcpy(&dataout[begin+sizeof(int)], &point.MBR[j].second, sizeof(int));
			// begin += 2*sizeof(int);
		// }
		// for(int j = 0; j<maxCap; j++){
			// memcpy(&dataout[begin], &point.children[j].id, sizeof(int));
			// begin += sizeof(int);
			// for(int k = 0; k<dim; k++){
				// memcpy(&dataout[begin], &point.children[j].MBR[k].first, sizeof(int));
				// memcpy(&dataout[begin+sizeof(int)], &point.children[j].MBR[k].second, sizeof(int));
				// begin += 2*sizeof(int);
			// }
			// memcpy(&dataout[begin+sizeof(int)], &point.children[j].MBR, 2*dim*sizeof(int));
		// }
		


		fhout->MarkDirty(point[0]);
		fhout->UnpinPage(point[0]);
		fhout->FlushPage(point[0]);
		foo(0, fhout);
		start += jump_size*dim;
	}
	foo(0, fhout);
	// cout << "Starting with Recursion\n";
	// recurse_until_root(fhout, start_pg_num, end_pg_num);	
	return ;
}

void insert(vector <int> *points, FileHandler *fhout){
	// To Implement
	return ;
}

// bool is_contained(vector< pair<int,int> > *mbr, vector <int> *points){
// 	for(int i = 0; i<mbr->size(); i++){
// 		if((points->at(i) > mbr->at(i).second) || (points->at(i) < mbr->at(i).first)){
// 			return false;
// 		}
// 	}
// 	return true;
// }

// bool is_leaf(Node *node){
// 	for(int i = 0; i<node->children.size(); i++){
// 		if(node->children[i].id != -1){
// 			return false; 
// 		}
// 	}
// 	return true;
// }

bool query(vector <int> *points, FileHandler *fhout, queue<int>*bfs_queue){
// 	PageHandler ph;
// 	// char* data;
// 	// Node *node;
// 	int pg_num;

// 	print_point(points);

// 	while(!bfs_queue->empty()){
// 		pg_num = bfs_queue->front();
// 		bfs_queue->pop();
// 		cout << "Node ID: " << pg_num << endl;
// 		ph = fhout->PageAt(pg_num);
// 		char *data = ph.GetData();
// 		Node node = Node();
// 		memcpy(&node, data, sizeof(Node));
// 		// (Node*) data;
// 		print_mbr(&node.MBR);
// 		if(!is_leaf(&node)){
// 			cout << "Internal\n";
// 			for(int i = 0; i<node.children.size(); i++){
// 				print_mbr(&node.children[i].MBR);
// 				if(is_contained(&node.children[i].MBR, points)){
// 					cout << "Something new in BFS\n";
// 					bfs_queue->push(node.children[i].id);
// 				}
// 			}
// 		}else{
// 			cout << "Leaf\n";
// 			if(is_contained(&node.MBR, points)){
// 				cout << "Found\n";
// 				return true;
// 			}
// 		}
// 		fhout->UnpinPage(pg_num);
// 	}
	return false;
}

int main(int argc, char *argv[]) {
	string inputfile = argv[1];
	maxCap = stoi(argv[2]);
	dim = stoi(argv[3]);
	string outputfile = argv[4];
	cout << inputfile << " " << maxCap << " " << dim << " " << outputfile << endl;

	nodesize = 1 + 2*dim + maxCap*(1+2*dim);
	// global_tree = Node(dim, maxCap);

	ifstream infile(inputfile);
	ofstream outfile (outputfile);
	FileHandler fhout = fm.CreateFile("bulkload.txt");

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
				insert(&points, &fhout);
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