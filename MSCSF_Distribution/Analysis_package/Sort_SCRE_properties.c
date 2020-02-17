#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[]){

    int i, j;
    FILE *in;
    FILE *out;
    int number_sims;
    
    double SR_value;
    double Cai_value;

    number_sims = atoi(argv[1]); 
    
    int spark[number_sims];
    float ti[number_sims];
    float tp[number_sims];
    float tf[number_sims];
    float RyR[number_sims];

    char *str_in, *str_out;
    str_in = malloc (50*sizeof(char));
    str_out = malloc (50*sizeof(char));
    str_in = argv[2];
    in = fopen(str_in, "r");

    sprintf(str_out, "Sorted_%s", str_in);
    out = fopen(str_out, "wt");

    // Read in original arrays
    for (i = 0; i < number_sims; i++)
    {
        fscanf(in, "%d %f %f %f %f\n", &spark[i], &ti[i], &tp[i], &tf[i], &RyR[i]); 
    } 

    // bubble sort by ti
    int tempi;
    float temp;
    
    for (j = 0; j < number_sims; j++)
    {

        for (i = 0; i < number_sims-1-j; i++)
        {
            if (ti[i] > ti[i+1])
            {
                tempi = spark[i];
                spark[i] = spark[i+1];
                spark[i+1] = tempi;

                temp = ti[i];
                ti[i] = ti[i+1];
                ti[i+1] = temp;

                temp = tp[i];
                tp[i] = tp[i+1];
                tp[i+1] = temp;

                temp = tf[i];
                tf[i] = tf[i+1];
                tf[i+1] = temp;

                temp = RyR[i];
                RyR[i] = RyR[i+1];
                RyR[i+1] = temp;
            }
        }  
    }

    for (i = 0; i < number_sims; i++)
    {   
        fprintf(out, "%d %f %f %f %f\n", spark[i], ti[i], tp[i], tf[i], RyR[i]); 
    }
}
