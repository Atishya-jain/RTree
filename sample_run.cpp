//Sample file for students to get their code running

#include<iostream>
#include<vector>
#include "file_manager.h"
#include "errors.h"
#include<cstring>

using namespace std;

struct temp{
	vector<int> b;
};

struct abc {
	int a = 10;
	temp c = temp();
};

int main() {
	FileManager fm;

	// Create a brand new file
	FileHandler fh = fm.CreateFile("temp.txt");
	cout << "File created " << endl;

	// Create a new page
	PageHandler ph = fh.NewPage ();
	char *data = ph.GetData ();

	abc d = abc();
	d.a = 100;
	d.c.b.push_back(200);

	// Store an integer at the very first location
	// int num = 5;
	memcpy (&data[0], &d, sizeof(abc));

	// // Store an integer at the second location
	// num = 1000;
	// memcpy (&data[4], &num, sizeof(int));

	// Flush the page
	fh.FlushPages ();
	cout << "Data written and flushed" << endl;

	// Close the file
	fm.CloseFile(fh);

	// Reopen the same file, but for reading this time
	fh = fm.OpenFile ("temp.txt");
	cout << "File opened" << endl;

	// Get the very first page and its data
	ph = fh.FirstPage ();
	data = ph.GetData ();

	// Output the first integer
	abc *e = (abc*) data;
	// memcpy (&e[0], &data[0], sizeof(abc));
	cout << "First number: " << e->a << " " << e->c.b[0] << endl;

	// // Output the second integer
	// memcpy (&num, &data[4], sizeof(int));
	// cout << "Second number: " << num << endl;;

	// Close the file and destory it
	fm.CloseFile (fh);
	fm.DestroyFile ("temp.txt");
}
