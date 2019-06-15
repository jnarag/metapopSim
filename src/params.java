import cern.jet.random.Exponential;
import cern.jet.random.Gamma;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

/**
 * Created by jayna on 16/05/2018.
 */
public class params {

    static int runTime;
    static int Npatches; // number of patches
    //static int I_total; // total number of infected cells in the simulation
    static int I; //
    static int S;
    static int R;
    static double gamma;
    static double beta; // per patch infectivity rate (1-20 per day in Lythgoe 2016)
    static double mig_blood_patch; // migration rate of viruses from blood to patch
    static double mig_patch_blood; // migration rate of virus from patch to blood
    static double U; // neutral mutation rate
    static double decay_in_blood; // decay rate of viruses in blood
    static double nu; // patch extinction rate
    static double c; // patch colonization rate


    static int seed;
    static RandomEngine randomGenerator;

    static Exponential extinctExpDist;
    static Exponential colonizeExpDist;


    public params() {

        runTime = 500;
        Npatches = 4;
        I = 5; // infected cells
        S = 25; // susceptible cells
        R = 0; // recovered cells
        beta = 0.2; // transmission rate
        gamma = 0.0; // recovery rate
        mig_blood_patch = 0.1 ; // migration rate
        //mig_patch_blood = 0.003; // migration rate
        //decay_in_blood = 0.2;
        nu = 0.05; //patch extinction rate;
        c = 0.02; // patch colonization rate;

        U = 0.1; // mutation rate
        seed = (int)Math.ceil(Math.random()*1000000000);
        randomGenerator = new MersenneTwister(1000);

        extinctExpDist = new Exponential(0.05, randomGenerator);
        colonizeExpDist = new Exponential(0.1, randomGenerator);


    }

    public int getNpatches() {

        return Npatches;
    }


}
