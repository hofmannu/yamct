/*
	main program only creating our software interface and starting it up
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 27.12.2020

	This looks quite dead for a main file. But what looks dead here is simply
	the beauty of object oriented programming.
*/

#include "interface.cuh"

int main(int *argcp, char**argv)
{
	interface GUI;
	GUI.InitWindow(argcp, argv);
	return 0;
}