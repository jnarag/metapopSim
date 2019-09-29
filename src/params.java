import cern.jet.random.Exponential;
import cern.jet.random.Gamma;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * Created by jayna on 16/05/2018.
 */
public class params {

    static double runTime;
    static double startTime;
    static int Npatches; // number of patches
    //static int I_total; // total number of infected cells in the simulation
    static int I; //
    static int S;
    static int R;
    static double gamma;
    static double beta; // per patch infectivity rate (1-20 per day in Lythgoe 2016)
    static double U; // neutral mutation rate
    static double nu; // patch extinction rate
    static double c; // patch colonization rate

    static int n_samples_per_time;


    static int seed;
    static RandomEngine randomGenerator;


    public params() {

        runTime = 5000;
        startTime = 0;
        Npatches = 3;
        I = 5; // infected cells
        S = 45; // susceptible cells
        R = 0; // recovered cells
        nu = 0.5; //patch extinction rate;
        c = 0.01; // patch colonization rate;

        beta = 1.0; // transmission rate
        gamma = 0.0; // recovery rate

        U = 10; // mutation rate
        seed = (int)Math.ceil(Math.random()*1000000000);
        randomGenerator = new MersenneTwister(seed);
        n_samples_per_time = 50;
    }

    public void print() {

        System.out.println("runTime "+params.runTime);
        System.out.println("extinction rate "+params.nu);
        System.out.println("colonization rate "+params.c);

    }


}
