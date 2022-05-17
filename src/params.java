import cern.jet.random.Exponential;
import cern.jet.random.Gamma;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * Created by jayna on 16/05/2018.
 */
public class params {


    static int n_sims;
    static double tau;
    static double runTime;
    static double startTime;
    static int Npatches; // number of patches
    static int I; //
    static int S;
    static String whModel; // 0 - logistic; 1 - SI model with death; 2 - SI model w/o death;
    static double beta; // infectivity rate (1-20 per day in Lythgoe 2016)
    static double death; // infected death rate (1-20 per day in Lythgoe 2016)
    static double r; //
    static double mu; // mutation rate per cell infection
    static Exponential benDFE; // beneficial fitness effect;
    static Exponential delDFE; // deleterious fitness effect;
    static double s_b;
    static double s_d;

    static double nu; // patch extinction rate
    static double c; // patch colonization rate

    static int n_samples_per_time;
    static int interval;


    static int seed;
    static RandomEngine randomGenerator;


    public params() {

        n_sims = 1;
        tau = 0.5;
        runTime = 5000;
        startTime = 1000;
        Npatches = 3;
        I = 5; // infected cells
        S = 45; // susceptible cells
        nu = 0.5; //patch extinction rate;
        c = 0.01; // patch colonization rate;

        beta = 0.1; // transmission rate;
        death = 0.0; // death rate of infected cells;
        r = 0.1; // within-patch growth rate;

        seed = (int)Math.ceil(Math.random()*1000);
        randomGenerator = new MersenneTwister(seed);
        n_samples_per_time = 10;
        interval = 10;
        whModel = "BSI";

        //convert mutation rate into per envelope gene per day (how many replications per day, how many errors per replication cycle)
        mu = 0.005; //0.005 = ~0.00005*2000*0.3/3 rate at which beneficial mutations appear at amino acid changing sites in envelope


        s_b = 0.01;
        s_d = 0.01;
        benDFE = new Exponential(1./s_b, randomGenerator);
        delDFE = new Exponential(1./s_d, randomGenerator);
    }

    public void print() {

        System.out.println("runTime "+params.runTime);
        System.out.println("N patches "+params.Npatches);
        System.out.println("patch size "+(params.S+params.I));
        System.out.println("extinction rate "+params.nu);
        //System.out.println("colonization rate "+params.c);
        System.out.println("beta "+params.beta);
        System.out.println("whModel "+params.whModel);
        System.out.println("mu " +params.mu);
        System.out.println("s_b " +params.s_b);
        System.out.println("s_d " +params.s_d);

        System.out.println("n_sims "+params.n_sims);



    }


}
