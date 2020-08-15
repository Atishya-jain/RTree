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
		for (int i = 0; i < dim; ++i){
			MBR.push_back(make_pair(INT_MIN, INT_MIN));
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
			MBR.push_back(make_pair(INT_MIN, INT_MIN));
		}
		for (int i = 0; i < maxCap; ++i){
			Child child = Child(dim);
			children.push_back(child);
		}
	}

	Node(int dim, int maxCap) : id(-1), pid(-1)
	{
		for (int i = 0; i < dim; ++i){
			MBR.push_back(make_pair(INT_MIN, INT_MIN));
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

void bulkload(string inputfile, int num_points){
	// Open the input file
	FileHandler fh = fm.OpenFile(inputfile.c_str());
	// Get first page
	PageHandler ph = fh.FirstPage();
	char *data = ph.GetData ();
	
	// Start your reading
	int id = 0;
	vector< Node > leaves;
	int start = 0;
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
		Node point = Node(id, -1, maxCap, dim);
		getMBRfrompoint(temp, &point.MBR);
		leaves.push_back(point);
		start += jump_size*dim;
		// cout << endl;
	}


	// To Implement
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

	if (infile.is_open() && outfile.is_open()){
		string line;
		while(getline(infile,line)){
			vector <string>commands;
			split_line(line, &commands);
			// cout << commands[0] << " " << commands.size() << '\n';
			if (commands[0] == "BULKLOAD"){
				bulkload(commands[1], stoi(commands[2]));
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