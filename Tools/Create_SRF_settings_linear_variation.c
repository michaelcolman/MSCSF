#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{

    double X, Y; // input value for distributions

    // timing
    double ti_sep, CF_ti_sep, k_ti_F1_ms, k_ti_F2_ms;
    double ti_A, ti, k_ti_F1, k_ti_F2;

    // duration/magnitude
    double MD, duration_width;

    double PSCR_threshold         = 1.0;
    double CaSR_max               = 2;
    double CaSR_min               = 1.0;
    double ti_sep_max             = 870;
    double ti_sep_min             = 30;
    double ti_width_max           = 500;
    double ti_width_min           = 20;
    double MD_max                 = 800;
    double MD_min                 = 50;
    double duration_width_max     = 250;
    double duration_width_min     = 30;
    double CaSR_width_H           = 1.0;

    double width_F1;

    char *string1 = (char*)malloc(500);
    FILE *out;

    for (Y = 0; Y < 1.05; Y +=0.05)
    {
        for (X = 0; X < 1.05; X +=0.05)
        {
            // ti based on X, duration based on Y
            // ti sep
            //ti_sep = (ti_sep_max - ti_sep_min)*exp(- (5*(X - CaSR_min)/(CaSR_max - CaSR_min)))  + ti_sep_min;
            //ti_sep = ti_sep_max - (ti_sep_max - ti_sep_min)*X;
            CF_ti_sep = 0.4;

            // ti widths
            //width_F1    =   (ti_width_max - ti_width_min)*pow((ti_sep - ti_sep_min)/(ti_sep_max - ti_sep_min), CaSR_width_H) + ti_width_min;
            width_F1    = ti_width_max - (ti_width_max - ti_width_min)*X;
            k_ti_F1    = width_F1;
            k_ti_F2    = 1.5 * k_ti_F1;

            ti_sep = 500+ti_width_max; // constant -> ensure that any ti is >500 ms, i.e. 500+max width = ti sep 

            // duration - MD
            //MD = (MD_max - MD_min)*exp(- (5*(Y - CaSR_min)/(CaSR_max - CaSR_min)))  + MD_min;
            //MD = MD_max - (MD_max - MD_min)*Y;
            //if (MD > MD_max) MD = MD_max;

            // duration width
            //duration_width = (duration_width_max - duration_width_min)*pow( (MD - MD_min)/(MD_max - MD_min)  , CaSR_width_H) + duration_width_min;
            duration_width = duration_width_max - (duration_width_max - duration_width_min)*Y;

            // make MD so that min duration of dist is always the same (30)
             MD = (150 - 125*Y) + duration_width;
            //MD = 50 + duration_width;
    
            //sprintf(string1, "SRF_settings_files_DC/SRF_settings_%0.2f_%0.2f_ti_%0.2f_%0.2f_D_%0.2f_%0.2f.txt", X, Y, ti_sep, k_ti_F1, MD, duration_width);
            sprintf(string1, "SRF_settings_files_DC/SRF_settings_%0.2f_%0.2f.txt", X, Y);
            out = fopen(string1, "wt");
            fprintf(out, "9\n");
            fprintf(out, "SRF_mode Direct_Control\n");
            fprintf(out, "SRF_Pset User_control\n");
            fprintf(out, "SRF_DC_PSCRE 1.000000\n");
            fprintf(out, "SRF_DC_CF %f\n", CF_ti_sep);
            fprintf(out, "SRF_DC_ti_sep %f\n", ti_sep);
            fprintf(out, "SRF_DC_ti_W1 %f\n", k_ti_F1);
            fprintf(out, "SRF_DC_ti_W2 %f\n", k_ti_F2);
            fprintf(out, "SRF_DC_MD %f\n", MD);
            fprintf(out, "SRF_DC_duration_W %f\n", duration_width);
            fclose(out);
        }
    }


}
