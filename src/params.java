import cern.jet.random.Exponential;
import cern.jet.random.Gamma;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * Created by jayna on 16/05/2018.
 */
public class params {

    static double tau;
    static double runTime;
    static double startTime;
    static int Npatches; // number of patches
    static int I; //
    static int S;
    static double beta; // infectivity rate (1-20 per day in Lythgoe 2016)
    static double death; // infected death rate (1-20 per day in Lythgoe 2016)
    static double r; // infected death rate (1-20 per day in Lythgoe 2016)

    static double nu; // patch extinction rate
    static double c; // patch colonization rate

    static int n_samples_per_time;
    static int interval;


    static int seed;
    static RandomEngine randomGenerator;

    static boolean het; //heterogeneity in patch growth rates;


    public params() {

        tau = 0.25;
        runTime = 5000;
        startTime = 0;
        Npatches = 3;
        I = 5; // infected cells
        S = 45; // susceptible cells
        nu = 0.5; //patch extinction rate;
        c = 0.01; // patch colonization rate;

        beta = 1.; // transmission rate;
        death = 1.; // death rate of infected cells;
        r = 0.1; // within-patch growth rate;

        seed = (int)Math.ceil(Math.random()*1000000000);
        randomGenerator = new MersenneTwister(seed);
        n_samples_per_time = 10;
        interval = 10;
        het = false;
    }

    public void print() {

        System.out.println("runTime "+params.runTime);
        System.out.println("N patches "+params.Npatches);
        System.out.println("patch size "+(params.S+params.I));
        System.out.println("extinction rate "+params.nu);
        System.out.println("colonization rate "+params.c);
        System.out.println("within patch growth rate "+params.r);
        System.out.println("patch heterogeneity "+params.het);


    }


}
