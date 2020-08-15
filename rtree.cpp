#include<iostream>
#include "file_manager.h"
#include "errors.h"
#include "constants.h"
#include<bits/stdc++.h>

using namespace std;

// Children inside a node
struct Child{
	int id;
	vector< pair<int,int>> MBR;	
	Child(int dim) : id(-1)
	{
		for (int i = 0; i < dim; i++){
			MBR.push_back(make_pair(INT_MAX, INT_MIN));
		}		
	}
};

// node of the tree
struct Node{
	int id;
	int pid;
	vector< pair<int,int> > MBR;
	vector< Child > children;
	Node(int myid, int parentid, int dim, int maxCap) : id(myid), pid(parentid)
	{
		for (int i = 0; i < dim; ++i){
			MBR.push_back(make_pair(INT_MAX, INT_MIN));
		}
		for (int i = 0; i < maxCap; ++i){
			Child child = Child(dim);
			children.push_back(child);
		}
	}

	Node(int dim, int maxCap) : id(-1), pid(-1)
	{
		for (int i = 0; i < dim; ++i){
			MBR.push_back(make_pair(INT_MAX, INT_MIN));
		}
		for (int i = 0; i < maxCap; ++i){
			Child child = Child(dim);
			children.push_back(child);
		}
	}
};

// Globals
// struct Node global_tree;
FileManager fm;
int jump_size = sizeof(int);
int dim, maxCap;
int root_page_number = INT_MIN;

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

void recurse_until_root(FileHandler *fh, int start, int end){
	cout << "Recursion: " << start << " " << end << endl;
	if(end > start){
		PageHandler ph;
		PageHandler ph_child;
		int start_pg_num = INT_MAX;
		int end_pg_num = INT_MIN;
		char* data;
		char* data_child;
		int counter = 0;
		while(true){
			ph = fh->NewPage();
			data = ph.GetData();
			Node point = Node(ph.GetPageNum(), -1, dim, maxCap);
			// getMBRChild(temp, &point.MBR);
			// Maintain start page as minimum page number 
			if(point.id <= start_pg_num){
				start_pg_num = point.id;
			}

			// Maintain end page as max page number 
			if(point.id >= end_pg_num){
				end_pg_num = point.id;
			}

			// Assign all children and set their parent id also
			counter = 0;
			while(start <= end){
				// Fetch this child page and copy into memory
				ph_child = fh->PageAt(start);
				data_child = ph_child.GetData();
				Node *child_point = (Node*) data_child;
				// memcpy(child_point, &data_child[0], sizeof(Node));
				// Setting parent's id and self child MBR
				child_point->pid = point.id;
				point.children[counter].MBR = child_point->MBR;
				point.children[counter].id = child_point->id;
				// Update your own MBR
				for(int i = 0; i<point.MBR.size(); i++){
					point.MBR[i].first = min(point.MBR[i].first, child_point->MBR[i].first);
					point.MBR[i].second = max(point.MBR[i].second, child_point->MBR[i].second);
				}
				// Copy back child_point into page and mark it dirty
				memcpy(&data_child[0], &child_point[0], sizeof(Node));
				fh->MarkDirty(child_point->id);
				fh->UnpinPage(child_point->id);
				counter++;
				start++;
				// If maxCap children done, break
				if(counter >= maxCap){
					break;
				}
			}
			memcpy(&data[0], &point, sizeof(Node));
			fh->MarkDirty(point.id);
			fh->UnpinPage(point.id);
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
	char *dataout;

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
		// cout << i << " ";
		vector<int> temp;
		for (int j = start; j<start+jump_size*dim; j+=4){
			int loc;
			memcpy(&loc, &data[j], jump_size);
			temp.push_back(loc);
			id += 1;
			// cout << loc << " ";
		}
		Node point = Node(id, -1, dim, maxCap);
		// Assign MBR as a point sized rectangle
		for(int i = 0; i<temp.size(); i++){
			point.MBR[i] = make_pair(temp[i], temp[i]);
		}

		// getMBRfrompoint(temp, &point.MBR);
		// leaves.push_back(point);
		phout = fhout->NewPage();
		dataout = phout.GetData();
		point.id = phout.GetPageNum();
		// Maintain start page as minimum page number 
		if(point.id <= start_pg_num){
			start_pg_num = point.id;
		}

		// Maintain end page as max page number 
		if(point.id >= end_pg_num){
			end_pg_num = point.id;
		}
		memcpy(&dataout[0], &point, sizeof(Node));
		// fh.FlushPage();
		fhout->MarkDirty(point.id);
		fhout->UnpinPage(point.id);
		start += jump_size*dim;
		// cout << endl;
	}
	cout << "Starting with Recursion\n";
	recurse_until_root(fhout, start_pg_num, end_pg_num);	
	return ;
}

void insert(vector <int> *points){
	// To Implement
	return ;
}

bool query(vector <int> *points){
	// To Implement
	return true;
}

int main(int argc, char *argv[]) {
	string inputfile = argv[1];
	maxCap = stoi(argv[2]);
	dim = stoi(argv[3]);
	string outputfile = argv[4];
	cout << inputfile << " " << maxCap << " " << dim << " " << outputfile << endl;

	// global_tree = Node(dim, maxCap);

	ifstream infile(inputfile);
	ofstream outfile (outputfile);
	FileHandler fhout = fm.CreateFile("bulkload.txt");

	if (infile.is_open() && outfile.is_open()){
		string line;
		while(getline(infile,line)){
			vector <string>commands;
			split_line(line, &commands);
			cout << commands[0] << " " << commands.size() << '\n';
			if (commands[0] == "BULKLOAD"){
				bulkload(commands[1], stoi(commands[2]), &fhout);
				outfile << "BULKLOAD\n\n\n";
			}else if(commands[0] == "INSERT"){
				vector<int> points;
				for (int i = 0; i < dim; i++){
					points.push_back(stoi(commands[i+1]));
				}
				insert(&points);
				outfile << "INSERT\n\n\n";
			}else{
				vector<int> points;
				for (int i = 0; i < dim; i++){
					points.push_back(stoi(commands[i+1]));
				}
				bool found = query(&points);				
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