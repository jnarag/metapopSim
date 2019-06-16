/**
 * Created by jayna on 17/05/2018.
 */

import cern.colt.list.DoubleArrayList;
import cern.jet.random.Normal;
import cern.jet.random.Poisson;
import cern.jet.random.Uniform;

import java.util.*;

public class onePatch_test {

    params params;
    infectionHistory patch_history;

    double tau; // time step
    int genotype;
    Map<Integer, List<Integer>> genotype_curr_per_patch;


    public onePatch_test() {

        params = new params();
        tau = 0.5;
        patch_history = new infectionHistory();

        genotype = 1;
        genotype_curr_per_patch = new HashMap<>();

    }


    public void runSim() {

        double t_curr = 0;
        double t_max = params.runTime;
        int Y_curr = params.I;
        int X_curr = params.S;
        int Z_curr = params.R;
        int total = Y_curr + X_curr + Z_curr;
        DoubleArrayList rates = new DoubleArrayList();
        List<Integer> curr_in_body = new ArrayList<>();
        List<Integer> persistent_genotypes = new ArrayList<>();


        int[] Y_curr_i = new int[params.Npatches];
        int[] X_curr_i = new int[params.Npatches];
        double[] timeThreshold_curr_i = new double[params.Npatches];
        double[] timeRefactory_curr_i = new double[params.Npatches];
        int[] threshold_i = new int[params.Npatches];


        //initializing patches with number of infected and susceptible cells

        Arrays.fill(Y_curr_i, 0);
        Arrays.fill(X_curr_i, (X_curr+Y_curr));
        Arrays.fill(timeThreshold_curr_i, 0.0);
        Arrays.fill(timeRefactory_curr_i, 0);

        int i1 = (int)Math.ceil(Uniform.staticNextDouble()*(double)Y_curr_i.length);
        int i2 = (int)Math.ceil(Uniform.staticNextDouble()*(double)Y_curr_i.length);

        Y_curr_i[i1-1] += Y_curr;
        X_curr_i[i1-1] -= Y_curr;
        Y_curr_i[i2-1] += Y_curr;
        X_curr_i[i2-1] -= Y_curr;


        initialiseHistory(patch_history, Y_curr_i, t_max);


        for (int p = 0; p < params.Npatches; p++) {

            List<Integer> genotype_curr_i = new ArrayList<>();

            for(Integer geno: patch_history.genotype) {

                int genotype_index = patch_history.genotype.indexOf(geno);

                if(patch_history.patch.get(genotype_index).equals(p)) {

                    int Y_temp = patch_history.prevalence.get(genotype_index).get(0);

                    int s = 0;
                    while(s < Y_temp) {
                        genotype_curr_i.add(geno);
                        s++;
                    }
                }
            }

            genotype_curr_per_patch.put(p, genotype_curr_i);

        }

        //controls the extinction rate...

        for(int i = 0; i < threshold_i.length; i++) {

            threshold_i[i] = Uniform.staticNextIntFromTo(20,30);

        }

        for (int t = 1; t < (t_max / tau); t++) {

            t_curr += tau;

            for (int i = 0; i < patch_history.prevalence.size(); i++) {

                List<Integer> prevalence = patch_history.prevalence.get(i);
                int prev_curr = prevalence.get(t - 1);
                patch_history.prevalence.get(i).set(t, prev_curr);

            }

            for (int i = 0; i < patch_history.prevalence.size(); i++) {

                Integer genotype_i = patch_history.genotype.get(i);

                if (patch_history.prevalence.get(i).get(t) > 0) {

                    Integer[] genotype_array = new Integer[patch_history.prevalence.get(i).get(t - 1)];
                    Arrays.fill(genotype_array, genotype_i);
                    curr_in_body.addAll(Arrays.asList(genotype_array));
                }

            }


            for (int p = 0; p < params.Npatches; p++) {

                int X_temp = X_curr_i[p];
                int Y_temp = Y_curr_i[p];


                rates.add(params.beta * X_temp * Y_temp / total); // infection rate
                rates.add(params.U * Y_temp); // mutation rate
                rates.add(((double) Y_temp / threshold_i[p]) * params.nu); // extinction
                rates.add((timeRefactory_curr_i[p]/t_curr) * params.c); // colonise

                Poisson poisson;
                int noOfRates = rates.size();
                int j;


                //events in patches that affect X and Y individuals
                for (int event = 0; event < noOfRates; event++) {

                    //curr_in_body.clear();

                    poisson = new Poisson(tau * rates.get(event), params.randomGenerator);
                    int num = poisson.nextInt();

                    int patch_index = p;

                    List<Integer> genotype_curr = new ArrayList<>();
                    genotype_curr.addAll(genotype_curr_per_patch.get(patch_index));


                    if (event == 0) { // new infection


                        int min = Math.min(X_curr_i[patch_index], num);

                        j = 0;

                        if (X_curr_i[patch_index] > 0 && Y_curr_i[patch_index] > 0) {

                            while (j < min) {

                                int r = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                                Integer chosen_genotype = genotype_curr.get(r);
                                int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

                                List<Integer> prevalence = new ArrayList<>();
                                prevalence.addAll(patch_history.prevalence.get(genotype_index));

                                int prevalence_curr = prevalence.get(t - 1);
                                int prevalence_next = prevalence.get(t);

                                X_curr_i[patch_index] = X_curr_i[patch_index] - 1;
                                Y_curr_i[patch_index] = Y_curr_i[patch_index] + 1;

                                if (Y_curr_i[patch_index] >= threshold_i[patch_index]) {

                                    timeThreshold_curr_i[patch_index] = t_curr;

                                }

                                if (prevalence_next > 0) {

                                    prevalence_curr = prevalence_next;
                                }

                                genotype_curr.add(chosen_genotype);
                                patch_history.prevalence.get(genotype_index).set(t, prevalence_curr + 1);

                                j++;
                            }

                        }

                    } else if (event == 1) { // mutations in patch_i

                        j = 0;

                        while (j < num && Y_curr_i[patch_index] > 0) {

                            //choose genotype
                            int index = (int) Math.floor(Uniform.staticNextDouble() * genotype_curr.size());
                            Integer chosen_genotype = genotype_curr.get(index);
                            int parent_index = patch_history.genotype.indexOf(chosen_genotype);


                            patch_history.logParent(chosen_genotype); // do we want to log parent index or genotype?
                            patch_history.logGenotype(genotype);
                            patch_history.logBirth(t_curr);
                            patch_history.logPatch(patch_index);
                            patch_history.logParentOrigin(0);

                            List<Integer> prevalence = new ArrayList<>();
                            prevalence.addAll(patch_history.prevalence.get(parent_index));

                            //prevalence is bit off need to check this
                            int prevalence_curr = prevalence.get(t - 1);
                            int prevalence_next = prevalence.get(t);

                            if (prevalence_next > 0) {

                                prevalence_curr = prevalence_next;
                            }

                            patch_history.prevalence.get(parent_index).set(t, prevalence_curr - 1);
                            List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[prevalence.size()]));
                            Collections.fill(new_prevalence, 0);
                            new_prevalence.set(t, 1);
                            patch_history.prevalence.add(new_prevalence);

                            genotype_curr.remove(chosen_genotype);
                            genotype_curr.add(genotype);
                            genotype++;

                            j++;
                        }


                    } else if (event == 2) {

                        // patch extinction

                        if (Y_curr_i[patch_index] >= threshold_i[patch_index]) {

                            Y_curr_i[patch_index] = 0;
                            X_curr_i[patch_index] = params.S + params.I;

                            //randomly choose a genotype to seed/colonize new patch

                            timeThreshold_curr_i[patch_index] = 0.0;
                            timeRefactory_curr_i[patch_index] = t_curr;

                            // track infection history (genotype_curr should be set to 0)
                            for (int i = 0; i < patch_history.prevalence.size(); i++) {

                                if (patch_history.patch.get(i) == patch_index) {
                                    patch_history.prevalence.get(i).set(t, 0);
                                }

                            }

                            threshold_i[patch_index] = Uniform.staticNextIntFromTo(20, 25);
                            //persistent_genotypes.addAll(chosen_genotypes); // only choose one genotype to persist from patch
                            // (needs to be from all patches)
                            genotype_curr.clear();

                        }
                    }
                    else {

                        // colonization

                        double p_colonize = 0;

                        if(curr_in_body.size() == 0) {


                        }

                        Collections.shuffle(curr_in_body);

                        if (Y_curr_i[patch_index] == 0) {

//                            p_colonize = params.extinctExpDist.cdf((t_curr - timeRefactory_curr_i[patch_index]));
                            double rand = Uniform.staticNextDouble();
//
                            if (rand < 0.5) {

                                // colonize

                                timeRefactory_curr_i[patch_index] = 0.0;

                                int newInfections = Uniform.staticNextIntFromTo(1, 6);

                                Y_curr_i[patch_index] = newInfections;
                                X_curr_i[patch_index] -= newInfections;

                                int r = (int) Uniform.staticNextDouble() * curr_in_body.size();
                                Integer colonizing_genotype = curr_in_body.get(r);

                                // need to go a few days back to see what genotypes could be lurking around.

                                for (int i = 0; i < newInfections; i++) {

                                    genotype_curr.add(genotype);

                                }
                                patch_history.logBirth(t_curr);
                                patch_history.logGenotype(genotype);
                                patch_history.logParent(colonizing_genotype);
                                patch_history.logPatch(patch_index);

                                List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                                Collections.fill(new_prevalence, 0);
                                new_prevalence.set(t, newInfections);
                                patch_history.prevalence.add(new_prevalence);
                                genotype++;
                                //curr_in_body.clear();

                            }

                        }

                    }

                    genotype_curr_per_patch.replace(patch_index, genotype_curr);

                }

                rates.removeAll(rates);

                //for (int i = 0; i < params.Npatches; i++) {
                    int size = 0;
                    for (int pr = 0; pr < patch_history.prevalence.size(); pr++) {

                        if (patch_history.patch.get(pr) == p) {
                            size += patch_history.prevalence.get(pr).get(t);
                        }
                    }

                    System.out.println("t_curr: " + t_curr +
                            ", toc_curr: " + timeThreshold_curr_i[p] +
                            ", tor_curr: " + timeRefactory_curr_i[p] +
                            ", threshold: " + threshold_i[p] +
                            ", Y_curr_i: " + Y_curr_i[p] +
                            ", X_curr_i: " + X_curr_i[p] +
                            ", genotype: " + (genotype - 1) +
                            ", totalInfections: " + size +
                            ", patch: " + p);


                //}

            }
            curr_in_body.clear();
            System.out.println();

        }
    }

    private int sumArray(int[] array) {

        int i = 0;
        int sum = 0;
        while(i < array.length) {

            sum += array[i];
            i++;
        }
        return sum;
    }

    private  Comparator<Integer>  comparator1 = new Comparator<Integer>() {
        @Override
        public int compare(Integer o1, Integer o2) {


            List<Integer> prevalence1 = patch_history.prevalence.get(o1);
            List<Integer> prevalence2 = patch_history.prevalence.get(o2);

            double diff = patch_history.getFitness(o1) - patch_history.getFitness(o2);
            if(diff == 0) {

                return 0;
            }
            else if(diff > 0) {
                return 1;
            }
            else{
                return -1;
            }
        }
    };

    public int chooseBin(double [] prob, double randomNo) {

        int bin = 0;

        while(prob[bin] < randomNo) {

            bin++;

            if(bin==prob.length){
                break;
            }

        }

        return bin;

    }

    private double getRandomNo() {

        Uniform randomNo = new Uniform(0,1, params.randomGenerator);

        return randomNo.nextDouble();

    }

    private void initialiseHistory(infectionHistory history, int [] prevalence, double t_max) {

        int n_patches = params.Npatches;
        int i = 0;
        while(i < n_patches) {
            history.logGenotype(genotype);
            history.logBirth(Double.NEGATIVE_INFINITY);
            history.logParent((int) Double.NEGATIVE_INFINITY);
            List<Integer> prevalence_through_time = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
            Collections.fill(prevalence_through_time, 0);
            prevalence_through_time.set(0, prevalence[i]);
            history.prevalence.add(prevalence_through_time);
            history.logPatch(i);
            genotype++;
            i++;
        }
    }



    public static void main(String [] args) {

        onePatch_test test = new onePatch_test();
        test.runSim();

    }

}
