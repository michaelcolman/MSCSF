#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>

using namespace std;

int main(int argc, char *argv[])
{

	int i, j, k;
	FILE *in;
	FILE *out;

	char *string1 = (char*)malloc(500);
	char *string2 = (char*)malloc(500);

	char * Ref;
	Ref = argv[1];

	char * Results_Ref_in = (char*)malloc(500);
    Results_Ref_in = argv[2];
	
    char * Results_Ref = (char*)malloc(500);

	int total_time;// = 5000;
	int threshold_time; // time during excitation, to be ignored

	int Nruns;

	total_time      = atoi(argv[3]);
	threshold_time  = atoi(argv[4]);
	Nruns  	    	= atoi(argv[5]);

	float time;
	float Vm;
	float Vm_max = -80;
	int start_switch = 0;
	float temp;

	int TA[3], TA_sum, TA_sum_sum;

	sprintf(string2, "%s_%s_focal_analysis.txt", Ref, Results_Ref_in);
	out = fopen(string2, "a");

	TA_sum_sum = 0;

	for (k = 0; k < Nruns; k++)
	{
		cout << "run = " << k << " of " << Nruns << endl;
		TA[0] = TA[1] = TA[2] = 0;
		start_switch = 0;

		// input file
		sprintf(Results_Ref, "%s_%04d", Results_Ref_in, k);

		sprintf(string1, "Outputs_0Dtissue_%s/Results_%s/CRU_cell1.txt", Ref, Results_Ref);


		in = fopen(string1, "r");
		if (!in)
		{
			cout << "File %s not found. skipping " << string1 << endl;
			TA[0] = -1;
		}
		else
		{
			printf("Reading filename: %s\n", string1);
			for (i = 0; i < total_time; i++)
			{
				for (j = 0; j < 25; j++)
				{
					if (j == 0) fscanf(in, "%f ", &time);
					else if (j == 1) fscanf(in, "%f ", &Vm);
					else if (j == 24) fscanf(in, "%f\n", &temp);
					else fscanf(in, "%f ", &temp);
				}

				if (time > threshold_time)
				{
					if (Vm > -20 && start_switch == 0)
					{
						start_switch = 1;
						TA[0] = 1;
					}

					if (Vm > Vm_max) Vm_max = Vm;

				}

			}
			fclose(in);
		} // end cell 1 else

		sprintf(string1, "Outputs_0Dtissue_%s/Results_%s/CRU_cell2.txt", Ref, Results_Ref);

		in = fopen(string1, "r");
		if (!in)
		{
			cout << "File not found. skipping " << endl;
			TA[1] = -1;
		}
		else
		{
			printf("Reading filename: %s\n", string1);
			for (i = 0; i < total_time; i++)
			{
				for (j = 0; j < 25; j++)
				{
					if (j == 0) fscanf(in, "%f ", &time);
					else if (j == 1) fscanf(in, "%f ", &Vm);
					else if (j == 24) fscanf(in, "%f\n", &temp);
					else fscanf(in, "%f ", &temp);
				}

				if (time > threshold_time)
				{
					if (Vm > -20 && start_switch == 0)
					{
						start_switch = 1;
						TA[1] = 1;
					}

					if (Vm > Vm_max) Vm_max = Vm;

				}

			}
			fclose(in);
		} // end cell 2 else


		sprintf(string1, "Outputs_0Dtissue_%s/Results_%s/CRU_cell3.txt", Ref, Results_Ref);

		in = fopen(string1, "r");
		if (!in)
		{
			cout << "File not found. skipping " << endl;
			TA[2] = -1;
		}
		else
		{
			printf("Reading filename: %s\n", string1);
			for (i = 0; i < total_time; i++)
			{
				for (j = 0; j < 25; j++)
				{
					if (j == 0) fscanf(in, "%f ", &time);
					else if (j == 1) fscanf(in, "%f ", &Vm);
					else if (j == 24) fscanf(in, "%f\n", &temp);
					else fscanf(in, "%f ", &temp);
				}

				if (time > threshold_time)
				{
					if (Vm > -20 && start_switch == 0)
					{
						start_switch = 1;
						TA[2] = 1;
					}

					if (Vm > Vm_max) Vm_max = Vm;

				}

			}
			fclose(in);
		} // end cell 3 else

		if (TA[0] > 0 || TA[1] > 0 || TA[2] > 0) TA_sum = 1;
		else TA_sum = 0;
		TA_sum_sum += TA_sum;

	} // end for

	fprintf(out, "%d %d %d\n", CaSR, TA_sum_sum, Nruns);
	fclose(out);
}
