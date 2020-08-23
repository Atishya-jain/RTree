//Sample file for students to get their code running

#include<iostream>
#include<vector>
#include "file_manager.h"
#include "errors.h"
#include<cstring>

using namespace std;

int main() {
	FileManager fm;

	// Create a brand new file
	FileHandler fh = fm.CreateFile("temp.txt");
	cout << "File created " << endl;
	int dim = 10;
	int points = 1000000;
	// Create a new page
	PageHandler ph;
	int number = 0;
	int start[dim];
	for (int i = 0; i < dim; ++i)
	{
		start[i] = 100;
	}
	while(true){
		ph = fh.NewPage();
		char *data = ph.GetData ();
		int pgnum = ph.GetPageNum();
		int offset = 0;
		while(offset + dim*4 < PAGE_CONTENT_SIZE){
			for(int i = 0; i<dim; i++){
				memcpy(&data[offset], &start[i], sizeof(int));
				start[i] += 100;
				offset+=4;
			}
			number++;
		}
		fh.MarkDirty(pgnum);
		fh.UnpinPage(pgnum);
		fh.FlushPage(pgnum);
		if(number >= points){
			break;
		}
	}
	for (int i = 0; i < dim; ++i)
	{
		cout << start[i] << ' ';
	}cout << endl;

	// // Get page number and say you just want to query something
	// int number = ph.GetPageNum();
	// fm.PrintBuffer();

	// // Unpin Page 
	// fh.UnpinPage(number);

	// Print Buffer (Buffer shows page unpinned correctly)
	// fm.PrintBuffer();

	// fh.FlushPage(number);

	// // Fetch page again and get data
	// ph = fh.PageAt(number);
	// data = ph.GetData();

	// Buffer still shows page unpinned ??
	// fm.PrintBuffer();
	// fh.UnpinPage(number);

	// ph = fh.PageAt(number);
	// data = ph.GetData();

	// // fh.UnpinPage(number);

	// fm.PrintBuffer();

	// // Flush the page
	// fh.FlushPages ();
	// cout << "Data written and flushed" << endl;

	// // Close the file
	// fm.CloseFile(fh);

	// // Reopen the same file, but for reading this time
	// fh = fm.OpenFile ("temp.txt");
	// cout << "File opened" << endl;

	// // Get the very first page and its data
	// ph = fh.FirstPage ();
	// data = ph.GetData ();

	// // Output the first integer
	// abc *e = (abc*) data;
	// // memcpy (&e[0], &data[0], sizeof(abc));
	// cout << "First number: " << e->a << " " << e->c.b[0] << endl;

	// // // Output the second integer
	// // memcpy (&num, &data[4], sizeof(int));
	// // cout << "Second number: " << num << endl;;

	// // Close the file and destory it
	fm.CloseFile (fh);
	// fm.DestroyFile ("temp.txt");
}
