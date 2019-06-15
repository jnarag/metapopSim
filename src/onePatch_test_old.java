/**
 * Created by jayna on 17/05/2018.
 */

import cern.colt.list.DoubleArrayList;
import cern.jet.random.Normal;
import cern.jet.random.Poisson;
import cern.jet.random.Uniform;

import java.util.*;

public class onePatch_test_old {

    params params;
    infectionHistory patch_history;

    double tau; // time step
    int genotype;
    Map<Integer, List<Integer>> genotype_curr_per_patch;


    public onePatch_test_old() {

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
        List<Integer> curr_in_blood = new ArrayList<>();
        List<Integer> persistent_genotypes = new ArrayList<>();


        int[] Y_curr_i = new int[params.Npatches];
        int[] X_curr_i = new int[params.Npatches];
        int[] Y_cum_i = new int[params.Npatches];
        double[] timeThreshold_curr_i = new double[params.Npatches];
        double[] timeRefactory_curr_i = new double[params.Npatches];


        //initializing patches with number of infected and susceptible cells

        Arrays.fill(Y_curr_i, Y_curr);
        Arrays.fill(X_curr_i, X_curr);
        Arrays.fill(Y_cum_i, Y_curr);
        Arrays.fill(timeThreshold_curr_i, Double.NaN);
        Arrays.fill(timeRefactory_curr_i, Double.NaN);


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

        int threshold = (int)Math.floor(Normal.staticNextDouble(params.S+params.I -10, 2));

        for (int t = 1; t <= (t_max / tau); t++) {


            // this is wrong, need to initialise genotype_curr per patch (genotype in diff patch cannot be a parent of new genotype or new infection)
            //Map<Integer, List<Integer>> genotype_curr_per_patch = new HashMap<>();
            //I think this has been dealt with 29/04/19

            t_curr+=tau;


            for (int i = 0; i < params.Npatches; i++) {

                if(Y_curr_i[i] == 0) {
                    System.out.println("x");
                }
                System.out.println("t_curr: " + t_curr + ", toc_curr: " + timeThreshold_curr_i[i] + ", tor_curr: " + timeRefactory_curr_i[i] +
                        ", threshold: " + threshold +
                        ", Y_curr_i: " + Y_curr_i[i] + ", X_curr_i: " + X_curr_i[i] + ", Y_cum_i: " + Y_cum_i[i] +
                        ", genotype: " + (genotype - 1) + ", patch: " + i + " curr_in_blood: " + curr_in_blood.size() + " genotype_curr: " + genotype_curr_per_patch.get(i));

                int X_temp = X_curr_i[i];
                int Y_temp = Y_curr_i[i];


                rates.add(params.beta * X_temp * Y_temp / total); // infection rate
                rates.add(params.U * Y_temp); // mutation rate
                rates.add(params.mig_blood_patch * Y_temp); // select from a pool of viruses from all patches
                rates.add(0.0);
                rates.add(0.0);
            }

            //System.out.println(rates);

            Poisson poisson;
            int noOfRates = rates.size();
            int j;


            //events in patches that affect X and Y individuals
            for (int event = 0; event < noOfRates; event++) {

                curr_in_blood.clear();

                poisson = new Poisson(tau * rates.get(event), params.randomGenerator);
                int num = poisson.nextInt();

                int patch_index = (int) Math.floor(event / noOfRates);

                List<Integer> genotype_curr = new ArrayList<>();
                genotype_curr.addAll(genotype_curr_per_patch.get(patch_index));


//                if(genotype_curr.size()==0) {
//
//                    System.out.println(genotype_curr);
//                    break;
//                }

                if (event % noOfRates == 0) { // new infection

                    int min = Math.min(X_curr_i[patch_index], num);
                    //Y_cum_i[patch_index] += min;

                    j = 0;

                    if (X_curr_i[patch_index] < min) {

                        int r = (int) Math.floor(getRandomNo() * genotype_curr.size());
                        int index = patch_history.genotype.indexOf(genotype_curr.get(r));

                        List<Integer> prevalence = new ArrayList<>();
                        prevalence.addAll(patch_history.prevalence.get(index));


                        int prevalence_curr = prevalence.get(t - 1);
                        int prevalence_next = prevalence.get(t);

                        if (prevalence_next > 0) {

                            prevalence_curr = prevalence_next;
                        }

                        patch_history.prevalence.get(index).set(t, prevalence_curr);


                    } else {
                        Y_cum_i[patch_index] += min;
                        //System.out.println(Y_curr_i[0]+ " "+ Y_curr_i[1]+" "+ Y_curr_i[2]);

//                        System.out.println("A "+X_curr_i[patch_index] + ": "+
//                                (X_curr_i[patch_index]-min)+ ", threshold: "+
//                        threshold);
                        while (j < min) {

                            int r = (int) Math.floor(getRandomNo() * genotype_curr.size());
                            Integer chosen_genotype = genotype_curr.get(r);
                            int index = patch_history.genotype.indexOf(chosen_genotype);

                            List<Integer> prevalence = new ArrayList<>();
                            prevalence.addAll(patch_history.prevalence.get(index));

                            //int time_index = t+1;

                            int prevalence_curr = prevalence.get(t - 1);
                            int prevalence_next = prevalence.get(t);


                            X_curr_i[patch_index] = X_curr_i[patch_index] - 1;
                            Y_curr_i[patch_index] = Y_curr_i[patch_index] + 1;

                            //int threshold = (int)Math.floor(Normal.staticNextDouble(params.S + params.I -5, 5));

                            if (Y_curr_i[patch_index] >= threshold) {

                                if(Double.isNaN(timeThreshold_curr_i[patch_index])) {

                                    timeThreshold_curr_i[patch_index] = t_curr;
                                }
                            }

                            if (prevalence_next > 0) {

                                prevalence_curr = prevalence_next;
                            }

                            genotype_curr.add(chosen_genotype);
                            patch_history.prevalence.get(index).set(t, prevalence_curr + 1);

                            //genotype_curr_per_patch.replace(patch_index, genotype_curr);

                            j++;
                        }

                    }

                }
//                else if (event % 4 == 1) {  // infected cell recovers
//
//                    int min = Math.min(Y_curr_i[patch_index], num);
//
//                    j = 0;
//
//                    if (Y_curr_i[patch_index] < min) {
//
//                        int r = (int) Math.floor(getRandomNo() * genotype_curr.size());
//                        int index = patch_history.genotype.indexOf(genotype_curr.get(r)); // index in the patch_history matrix
//                        List<Integer> prevalence = new ArrayList<>();
//                        prevalence.addAll(patch_history.prevalence.get(index));
//
//
//                        int prevalence_curr = prevalence.get(t - 1);
//                        int prevalence_next = prevalence.get(t);
//
//                        if(prevalence_next > 0) {
//
//                            prevalence_curr = prevalence_next;
//                        }
//                        patch_history.prevalence.get(index).set(t, prevalence_curr);
//
//                    } else {
//
//                        while (j < min) {
//                            int r = (int) Math.floor(getRandomNo() * genotype_curr.size());
//
//                            int index = patch_history.genotype.indexOf(genotype_curr.get(r));
//                            List<Integer> prevalence = new ArrayList<>();
//                            prevalence.addAll(patch_history.prevalence.get(index));
//
//                            genotype_curr.remove(r);
//
//
//                            Z_curr_i[patch_index] = Z_curr_i[patch_index] + 1;
//                            Y_curr_i[patch_index] = Y_curr_i[patch_index] - 1;
//
//
//                            int prevalence_curr = prevalence.get(t - 1);
//                            int prevalence_next = prevalence.get(t);
//
//                            if(prevalence_next > 0) {
//                                prevalence_curr = prevalence_next;
//                            }
//
//                            patch_history.prevalence.get(index).set(t, prevalence_curr-1);
//
////                            if(prevalence_curr != genotype_curr.size()) {
////                                System.out.println("bbb");
////                            }
//
//                            j++;
//                        }
//                    }
//                }
                else if (event % noOfRates == 1) { // mutations in patch_i

                    j = 0;

                    List<Integer> genotype_curr_copy = new ArrayList<>();
                    genotype_curr_copy.addAll(genotype_curr);

//                    System.out.println("1x " +genotype_curr_per_patch);
//                    System.out.println("1x " +Y_curr_i[patch_index]);
//
//                    System.out.println("1x " +", "+X_curr_i[patch_index]);
//                    System.out.println("1x " +", "+curr_in_blood);



                    while (j < num) {

                        int index = (int) Math.floor(getRandomNo() * genotype_curr.size());

                        int parent_index = patch_history.genotype.indexOf(genotype_curr.get(index));

                        patch_history.logParent(genotype_curr.get(index)); // do we want to log parent index or genotype?
                        patch_history.logGenotype(genotype);
                        patch_history.logBirth(t_curr);
                        patch_history.logPatch(patch_index);
                        patch_history.logParentOrigin(0);

                        List<Integer> prevalence = new ArrayList<>();
                        prevalence.addAll(patch_history.prevalence.get(parent_index));

                        genotype_curr_copy.remove(genotype_curr.get(index));

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

                        genotype_curr_copy.add(genotype);
                        genotype++;

                        //genotype_curr_per_patch.replace(patch_index, genotype_curr);
                        j++;
                    }

                    genotype_curr.clear();
                    genotype_curr.addAll(genotype_curr_copy);

                }
                else if (event % noOfRates == 2) {

                    int min = Math.min(X_curr_i[patch_index], num);

                    // viruses from blood to patch
                    if (min < X_curr_i[patch_index]) {

                        j = 0;

                        Y_cum_i[patch_index] += num;

                        while (j < min) {

                            Set<Integer> patchIndices = genotype_curr_per_patch.keySet();

                            for (Integer i : patchIndices) {

                                curr_in_blood.addAll(genotype_curr_per_patch.get(i));
                            }

                            if (Y_curr_i[patch_index] >= threshold) {

                                if(Double.isNaN(timeThreshold_curr_i[patch_index])) {

                                    timeThreshold_curr_i[patch_index] = t_curr;
                                }
                            }

                            Collections.sort(curr_in_blood);

                            int r = (int) Math.floor(getRandomNo() * curr_in_blood.size());

                            int parent = curr_in_blood.get(r);
                            genotype_curr.add(genotype);
                            patch_history.logGenotype(genotype);
                            patch_history.logParent(parent);
                            patch_history.logBirth(t_curr);
                            List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau) + 2]));
                            Collections.fill(new_prevalence, 0);
                            new_prevalence.set(t, 1);
                            patch_history.prevalence.add(new_prevalence);

                            X_curr_i[patch_index] = X_curr_i[patch_index] - 1;
                            Y_curr_i[patch_index] = Y_curr_i[patch_index] + 1;
                            j++;
                            genotype++;
                        }
                    }

                }

                else if (event % noOfRates == 3) {
                    // patch extinction

                    // does the patch go extinct?

                    double p_extinct = 0;
                    if (!Double.isNaN((timeThreshold_curr_i[patch_index]))) {
//
                        p_extinct = params.extinctExpDist.cdf((t_curr - timeThreshold_curr_i[patch_index]));
                    }
//
//                        //System.out.println(Math.random()+ ", " + (t_curr - timeThreshold_curr_i[patch_index]) + ", " + p_extinct);
//
                    if (Math.random() < p_extinct) {

                        // patch goes extinct

                        Y_curr_i[patch_index] = 0;
                        X_curr_i[patch_index] = params.S + params.I;
                        Y_cum_i[patch_index] = 0;

                        //randomly choose a genotype to seed/colonize new patch
                        int r = (int) Math.floor(getRandomNo() * genotype_curr.size());
                        Integer chosen_genotype = genotype_curr.get(r);

                        timeThreshold_curr_i[patch_index] = Double.NaN;
                        timeRefactory_curr_i[patch_index] = t_curr;

//
                        for (Integer g : genotype_curr) {


                            int g_index = patch_history.genotype.indexOf(g);

                            if (g_index == r) {


                                patch_history.prevalence.get(g_index).set(t, 1);

                            } else {
                                patch_history.prevalence.get(g_index).set(t, 0);
                            }
                        }


                        threshold = (int) Math.ceil(Normal.staticNextDouble((params.S+params.I-10) , 2));
                        persistent_genotypes.add(chosen_genotype);
                        genotype_curr.clear();
                        //genotype_curr.add(chosen_genotype);

                    }
                }
                else {

                    // colonization

                    double p_colonize = 0;
                    if (!Double.isNaN((timeRefactory_curr_i[patch_index]))) {

                        p_colonize = params.extinctExpDist.cdf((t_curr - timeRefactory_curr_i[patch_index]));

                        //System.out.println(">>"+p_colonize + ", "+ (t_curr - timeRefactory_curr_i[patch_index]));
                        if (Math.random() < p_colonize) {

                            // colonize

                            timeRefactory_curr_i[patch_index] = Double.NaN;

                            int newInfections = (int) Math.ceil(Math.exp(Normal.staticNextDouble(0.1,1)));
                            Y_curr_i[patch_index] = newInfections;
                            X_curr_i[patch_index] -= newInfections;
                            Y_cum_i[patch_index] = newInfections;

                            Integer colonizing_genotype = persistent_genotypes.get(0);

                            for (int i = 0; i < newInfections; i++) {

                                genotype_curr.add(genotype);


                            }
                            patch_history.logBirth(t_curr);
                            patch_history.logGenotype(genotype);
                            patch_history.logParent(colonizing_genotype);

                            List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau) + 2]));
                            Collections.fill(new_prevalence, 0);
                            new_prevalence.set(t, 1);
                            patch_history.prevalence.add(new_prevalence);
                            genotype++;

                        }
                    }
                }
//                }

//                System.out.println(genotype_curr);
//                System.out.println(genotype_curr_per_patch);

                genotype_curr_per_patch.replace(patch_index, genotype_curr);
//                System.out.println(genotype_curr_per_patch);

            }

            rates.removeAll(rates);


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
            List<Integer> prevalence_through_time = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau) + 2]));
            Collections.fill(prevalence_through_time, 0);
            prevalence_through_time.set(0, prevalence[i]);
            history.prevalence.add(prevalence_through_time);
            history.logPatch(i);
            genotype++;
            i++;
        }
    }



    public static void main(String [] args) {

        onePatch_test_old test = new onePatch_test_old();
        test.runSim();
    }

}
