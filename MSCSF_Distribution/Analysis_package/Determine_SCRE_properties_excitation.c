#include <stdio.h>

int main(int argc, char *argv[])
{

    int i, j;
    FILE *in;
    FILE *out;

    char * filename_in;
    filename_in = argv[1];

    in = fopen(filename_in, "r");
    printf("Reading filename: %s\n", filename_in);
    out = fopen("SCRE_properties_log.dat", "a");

    int total_time;// = 5000;
    int threshold_time; // time during excitation, to be ignored
    int BCL;
    int beats;

    BCL = atoi(argv[2]);
    beats = atoi(argv[3]);
    total_time = atoi(argv[4]);

    threshold_time = (beats-1)*BCL + 150;

    printf("BCL = %d\nbeats = %d\nTotal_time =%d\nThreshold time = %d\n", BCL, beats, total_time, threshold_time);

    float time;
    float RyRo;
    float RyRo_max = 0;
    float ti; // initiation time
    float tp; // peak time
    float tf; // finish time
    int spark_yes = 0;
    int start_switch = 0;
    int end_switch = 0;
    float temp;
    
    for (i = 0; i < total_time; i++)
    {
       for (j = 0; j < 25; j++)
       {
           if (j == 0) fscanf(in, "%f ", &time);
           else if (j == 8) fscanf(in, "%f ", &RyRo);
           else if (j == 24) fscanf(in, "%f\n", &temp);
           else fscanf(in, "%f ", &temp);
       }

        if (time > threshold_time)
        {
            if (RyRo > 0.01 && spark_yes == 0 && start_switch == 0) // calcs time of intial deflection, IF later turns into full wave
            {
                ti = time;
                start_switch = 1;
            }

            if (start_switch == 1 && RyRo < 0.01 && spark_yes == 0) start_switch = 0; // resets if deflection goes small again before full wave

            if (RyRo > 0.035 && spark_yes == 0 && start_switch == 1) // actual wave
            {
                spark_yes = 1;
            }

            if (spark_yes == 1 && end_switch == 0)
            {
                if (RyRo > RyRo_max) // peak time and RyR o max value
                {
                    RyRo_max = RyRo;
                    tp = time;
                }

                if (RyRo < 0.01) // end wave
                {
                    tf = time;
                    end_switch = 1;
                }
            }
        }

    }

    if (spark_yes != 1) { ti = tf = tp = RyRo_max = 0; }

    fprintf(out, "%d %f %f %f %f\n", spark_yes, ti, tp, tf, RyRo_max);
    fclose(out);

}
