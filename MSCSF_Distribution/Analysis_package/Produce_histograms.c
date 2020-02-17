#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[]){

    int i;
    FILE *in;
    FILE *out;
    FILE *out2;
    int number_sims;
    int total_time;
    int spark_yes;
    int number_SCRE = 0;
    
    number_sims = atoi(argv[2]); 
    total_time = atoi(argv[3]);
    
    int Hist[total_time/5];
    int Hist_duration[1500/5];
    int time_range;
    int duration;
    int ti_int;

    for (i = 0; i < total_time/5; i++) Hist[i] = 0;
    for (i = 0; i < 1500/5; i++) Hist_duration[i] = 0;

    char *str_in;
    str_in = malloc (50*sizeof(char));
    str_in = argv[1];
    in = fopen(str_in, "r");
    printf("Filename %s read\n", str_in);

    out2 = fopen("Summary_data.dat", "a");
    fprintf(out2, "%s ", str_in);

    float ti, tp, tf, RyRo;
    int first_spark_switch = 0;
    int tp_non_zero;

    for (i = 0; i < number_sims; i++)
    {
        fscanf(in, "%d %f %f %f %f\n", &spark_yes, &ti, &tp, &tf, &RyRo); 
        if (spark_yes == 1)
        {
            if (first_spark_switch == 0)
            {   
                tp_non_zero = ti;
                first_spark_switch = 1;
            }       
            number_SCRE ++;
            time_range = ti/5;
            Hist[time_range]++;
            duration = (tf-ti)/5;
            Hist_duration[duration]++;
        }
    } 

    printf("Time of earliest SCRE = %d\n", tp_non_zero);
    fprintf(out2, "%d ", tp_non_zero);

    fclose(in);
    fclose(out2);

    out = fopen("Initiation_time_histogram.dat", "wt");

    for (i = 0; i < total_time/5; i++)
    {
        fprintf(out, "%d %d\n", (i*5)+2, Hist[i]);
    }
    fclose(out);
    
    out = fopen("Duration_histogram.dat", "wt");

    for (i = 0; i < 1500/5; i++)
    {
        fprintf(out, "%d %d\n", (i*5)+2, Hist_duration[i]);
    }
    fclose(out);
}
