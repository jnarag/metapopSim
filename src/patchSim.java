/**
 * Created by jayna on 17/05/2018.
 */

import cern.colt.list.DoubleArrayList;
import cern.jet.random.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

import org.eclipse.collections.api.list.primitive.MutableIntList;
import org.eclipse.collections.impl.factory.primitive.IntLists;
import utils.*;


public class patchSim {

    params params;
    infectionHistory patch_history;

    int genotype;
    Map<Integer, List<Integer>> genotype_curr_per_patch;
    List<Integer> emptyPatchList;
    List<Integer> occupyPatchList;
    List<Integer> thresholdList;

    int total_infected;
    int total_genotypes;
    double threshold;
    double r_i;

    public patchSim() {

        patch_history = new infectionHistory();
        genotype = 1;
        genotype_curr_per_patch = new HashMap<>();
        emptyPatchList = new ArrayList<>();
        occupyPatchList = new ArrayList<>();
        thresholdList = new ArrayList<>();

        total_genotypes = 0;
        total_infected = 0;
        r_i = params.r;

    }


    public void runSim(params params) {


        double tau = params.tau;
        boolean writeOutput = false;

        String outputfile = "patchSim_extinction_" + params.nu + "_col_" + params.c +
                "_r_" + params.r + "_Npatches_" + params.Npatches +
                "_patchSize_" + (params.S + params.I + ".txt");


        String summaryFile = "summary_patchSim_extinction_" + params.nu + "_col_" + params.c +
                "_r_" + params.r + "_Npatches_" + params.Npatches +
                "_patchSize_" + (params.S + params.I + ".txt");

        FileWriter writer1 = null;
        FileWriter writer2 = null;
        if(writeOutput) {
            try {
                writer1 = new FileWriter(new File(outputfile));
                writer1.write("Time\tPatch_no\tTotal_infected\tUnique_genotypes\n");

                writer2 = new FileWriter(new File(summaryFile));
                writer2.write("Time\tTotal_infected\tUnique_genotypes\tOccupied_patches\tGenealogical_diversity\tTMRCA\n");
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

        double t_curr = tau;
        double t_max = params.runTime;
        int N = params.S + params.I;

        DoubleArrayList betweenPatchRates = new DoubleArrayList();


        double[] timeRefactory_curr_i = new double[params.Npatches];
        threshold = 0.5*N;


        //initializing patches with number of infected and susceptible cells

        Arrays.fill(timeRefactory_curr_i, -5);


        double colonize = params.c;
        double extinct = params.nu;

        initialiseHistory(patch_history, params.I, t_max);


        //IntStream.range(0, params.Npatches).forEach(i -> genotype_curr_per_patch.put(i, new ArrayList<>()));

        for(Integer genotype: patch_history.genotype) {
            int genotype_index = patch_history.genotype.indexOf(genotype);

            int patch = patch_history.getPatch(genotype_index);
            occupyPatchList.add(patch);

            int n_infected = patch_history.prevalence.get(genotype_index).get(0);
            List<Integer> Y_temp = new ArrayList<>( Arrays.asList(new Integer[n_infected]));

            Collections.fill(Y_temp, genotype);

            genotype_curr_per_patch.put(patch, Y_temp);

            total_infected += n_infected;
            total_genotypes += 1;

        }

        MutableIntList listNew =
                IntLists.mutable.withAll(
                        IntStream.rangeClosed(1, 100));



        String sim = "ext_"+ extinct+"_col_"+colonize+"_Npatch_"+ params.Npatches +"_PatchSize_"+(N);

        int N_occupied = occupyPatchList.size();


        for (int t = 1; t < (t_max / tau); t++) {


//            //determine what is present currently across all patches
//            update_curr_in_body(curr_in_body, t);





            if(t_curr%20==0) {

                updatePatchStatistics(writer1, t_curr, writeOutput);

                List<Double> globalDiversity = updateDiversity(false);
                double tmrca = globalDiversity.get(1);

                System.out.println("t_curr \t" + t_curr + "\ttotal_infected: \t" + Math.log10(total_infected) + "\ttotal_genotypes: \t" + total_genotypes +
                        "\tpatch_occupied: \t"+ N_occupied + "\t"+(tmrca));


                if(writeOutput) {


                    writeOutput(writer2, t_curr, total_infected, total_genotypes, N_occupied, globalDiversity);
                }

            }

            betweenPatchRates.add(extinct*thresholdList.size());
            betweenPatchRates.add(colonize*(params.Npatches-N_occupied));
            int noOfRates = betweenPatchRates.size();

            Poisson poisson;

            for(int event = 0; event < noOfRates; event++) {

                int j;

                poisson = new Poisson(betweenPatchRates.get(event)*tau, params.randomGenerator);
                int num = poisson.nextInt();

                if(num == 0) continue;

                //extinction

                switch(event) {
                    case 0:

                        int min = Math.min(thresholdList.size(), num);

                        j = 0;

                        while (j < min) {

                            // patch extinction

                            int rand = Uniform.staticNextIntFromTo(0, thresholdList.size() - 1);
                            int patch_index = thresholdList.get(rand);

                            // track infection history (genotype_curr should be set to 0)

                            Integer g = genotype_curr_per_patch.get(patch_index).get(0);


                            int index = patch_history.genotype.indexOf(g);

                            patch_history.prevalence.get(index).set(t, 0);
                            patch_history.death.set(index, t_curr);

                            total_infected -= genotype_curr_per_patch.get(patch_index).size();
                            total_genotypes -= 1;


                            genotype_curr_per_patch.remove((Integer)patch_index);
                            occupyPatchList.remove((Integer) patch_index);
                            thresholdList.remove((Integer) patch_index);
                            timeRefactory_curr_i[patch_index] = t_curr;

                            j++;

                        }
                        N_occupied = occupyPatchList.size();
                        break;

                        // patch colonisation


                    case 1:

                        min = Math.min((params.Npatches-N_occupied), num);

                        j = 0;

                        while (j < min) {

                            j++;
                            int patch_index = Uniform.staticNextIntFromTo(0, params.Npatches - 1);

                            if(occupyPatchList.contains(patch_index)) {
                                while(occupyPatchList.contains(patch_index)) {
                                    patch_index = Uniform.staticNextIntFromTo(0, params.Npatches - 1);
                                }
                            }
                            //int patch_index = emptyPatchList.get(0);

                            if(!genotype_curr_per_patch.containsKey(patch_index)) {

                                genotype_curr_per_patch.put(patch_index, new ArrayList<>());
                            }

                            if (total_infected == 0 || occupyPatchList.size() == 0) {
                                break;
                            }

                            if (t_curr - timeRefactory_curr_i[patch_index] < 5) {
                                continue;
                            }

                            int s = Uniform.staticNextIntFromTo(0, occupyPatchList.size() - 1);

                            int source = occupyPatchList.get(s);


                            int newInfections = Uniform.staticNextIntFromTo(1, 3);


                            Integer colonizing_genotype = genotype_curr_per_patch.get(source).get(0);  //curr_in_body.get(r);

                            for (int i = 0; i < newInfections; i++) {

                                genotype_curr_per_patch.get(patch_index).add(genotype);

                            }

                            total_infected += newInfections;
                            total_genotypes += 1;

                            patch_history.logBirth(t_curr);
                            patch_history.logDeath(Double.POSITIVE_INFINITY);
                            patch_history.logGenotype(genotype);
                            patch_history.logParent(colonizing_genotype);
                            patch_history.logPatch(patch_index);

                            List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                            Collections.fill(new_prevalence, 0);
                            new_prevalence.set(t, newInfections);
                            patch_history.prevalence.add(new_prevalence);
                            genotype++;

                            occupyPatchList.add(patch_index);

                        }
                        N_occupied = occupyPatchList.size();
                        break;
                }

            }

            betweenPatchRates.clear();


            for (Integer p: genotype_curr_per_patch.keySet()) {

                List<Integer> genotype_curr = new ArrayList<>();
                genotype_curr.addAll(genotype_curr_per_patch.get(p));

                int Y_temp = genotype_curr.size();
                if(Y_temp == 0) continue;

                double infection_rate = ((r_i*Y_temp)*(1.0-(Y_temp*1.0/N*1.0)));

                poisson = new Poisson(infection_rate*tau, params.randomGenerator);
                int num = poisson.nextInt();

                int min = Math.min((N - Y_temp), num); // max number of new infected cells is bounded by number of uninfected cells


                Integer chosen_genotype = genotype_curr.get(0);
                int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

                int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t - 1);
                int prevalence_next = patch_history.prevalence.get(genotype_index).get(t);


                if (prevalence_next > 0) {

                    prevalence_curr = prevalence_next;
                }

                List<Integer> temp_list = new ArrayList<>(Arrays.asList(new Integer[min]));
                Collections.fill(temp_list, chosen_genotype);
                genotype_curr.addAll(temp_list);
                patch_history.prevalence.get(genotype_index).set(t, prevalence_curr + min);

                if ((prevalence_curr + min) >= threshold && !thresholdList.contains(new Integer(p))) {

                    thresholdList.add(p);
                }

                total_infected += min;
                genotype_curr_per_patch.replace(p, genotype_curr);
            }

            t_curr += tau;

        }

        if(writeOutput) {
            try {
                writer1.close();
                writer2.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void runSimExpGrowth(params params) {

        double tau = params.tau;
        boolean writeOutput = false;

        String outputfile = "patchSim_extinction_" + params.nu + "_col_" + params.c +
                "_r_" + params.r + "_Npatches_" + params.Npatches +
                "_patchSize_" + (params.S + params.I + ".txt");


        String summaryFile = "summary_patchSim_extinction_" + params.nu + "_col_" + params.c +
                "_r_" + params.r + "_Npatches_" + params.Npatches +
                "_patchSize_" + (params.S + params.I + ".txt");

        FileWriter writer1 = null;
        FileWriter writer2 = null;
        if(writeOutput) {
            try {
                writer1 = new FileWriter(new File(outputfile));
                writer1.write("Time\tPatch_no\tTotal_infected\tUnique_genotypes\tGenealogical_diversity\tTMRCA\n");

                writer2 = new FileWriter(new File(summaryFile));
                writer2.write("Time\tTotal_infected\tUnique_genotypes\tOccupied_patches\tGenealogical_diversity\tTMRCA\n");
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

        double t_curr = tau;
        double t_max = params.runTime;
        int N = params.S + params.I;

        DoubleArrayList withinPatchRates = new DoubleArrayList();
        DoubleArrayList betweenPatchRates = new DoubleArrayList();
        List<Integer> curr_in_body = new ArrayList<>();


        //double[] beta_i = new double[params.Npatches];
        double[] timeRefactory_curr_i = new double[params.Npatches];
        threshold = 0.5*N;


        //initializing patches with number of infected and susceptible cells

        Arrays.fill(timeRefactory_curr_i, -5);

        double colonize = params.c;
        double extinct = params.nu;
        double r = params.r;


        initialiseHistory(patch_history, params.I, t_max);



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

        for(int i = 0; i < params.Npatches; i++) {

            //beta_i[i] = beta;

            if(genotype_curr_per_patch.get(i).size() == 0) {
                emptyPatchList.add(i);
            }
            else{
                occupyPatchList.add(i);
            }
        }

        String sim = "ext_"+ extinct+"_col_"+colonize+"_Npatch_"+ params.Npatches +"_PatchSize_"+(N);

        for (int t = 1; t < (t_max / tau); t++) {

            t_curr += tau;


            //determine what is present currently across all patches
            update_curr_in_body(curr_in_body, t);


            if(t_curr%20==0) {

                updatePatchStatistics(writer1, t_curr, writeOutput);

                List<Double> globalDiversity = updateDiversity(false);
                double tmrca = globalDiversity.get(1);

                System.out.println("t_curr \t" + t_curr + "\ttotal_infected: \t" + Math.log10(total_infected) + "\ttotal_genotypes: \t" + total_genotypes +
                        "\tpatch_occupied: \t"+ occupyPatchList.size() + "\t"+(tmrca));


                if(writeOutput) {
                    writeOutput(writer2, t_curr, total_infected, total_genotypes, occupyPatchList.size(), globalDiversity);
                }

            }

            betweenPatchRates.add(extinct*thresholdList.size());
            betweenPatchRates.add(colonize*emptyPatchList.size());
            int noOfRates = betweenPatchRates.size();

            Poisson poisson;

            for(int event = 0; event < noOfRates; event++) {

                int j;

                poisson = new Poisson(betweenPatchRates.get(event)*tau, params.randomGenerator);
                int num = poisson.nextInt();

                if(num == 0) continue;

                //extinction

                switch(event) {
                    case 0:

                        int min = Math.min(thresholdList.size(), num);

                        j = 0;
                        while (j < min) {

                            // patch extinction

                            int rand = Uniform.staticNextIntFromTo(0, thresholdList.size() - 1);
                            int patch_index = thresholdList.get(rand);

                            // track infection history (genotype_curr should be set to 0)

                            Set<Integer> uniqueGenotypes = new HashSet<>();
                            uniqueGenotypes.addAll(genotype_curr_per_patch.get(patch_index));

                            for (Integer g : uniqueGenotypes) {

                                int index = patch_history.genotype.indexOf(g);

                                if (patch_history.patch.get(index) == patch_index) {
                                    patch_history.prevalence.get(index).set(t, 0);
                                    patch_history.death.set(index, t_curr);

                                }

                            }

                            genotype_curr_per_patch.get(patch_index).clear();
                            occupyPatchList.remove((Integer) patch_index);
                            thresholdList.remove((Integer) patch_index);
                            emptyPatchList.add(patch_index);
                            timeRefactory_curr_i[patch_index] = t_curr;

                            j++;

                        }
                        break;



                        // patch colonisation


                    case 1:

                        min = Math.min(emptyPatchList.size(), num);

                        j = 0;
                        while (j < min) {

                            int rand = Uniform.staticNextIntFromTo(0, emptyPatchList.size() - 1);

                            int patch_index = emptyPatchList.get(rand);

                            if (t_curr - timeRefactory_curr_i[patch_index] < 5) {
                                j++;
                                continue;
                            }

                            if (curr_in_body.size() == 0 || occupyPatchList.size() == 0) {

                                break;
                            }

                            int s = Uniform.staticNextIntFromTo(0, occupyPatchList.size() - 1);
                            if (occupyPatchList.size() == 0) {
                                System.out.println();
                            }
                            int source = occupyPatchList.get(s);


                            int newInfections = Uniform.staticNextIntFromTo(1, 3);

                            //Y_curr_i[patch_index] = newInfections;

                            Integer colonizing_genotype = genotype_curr_per_patch.get(source).get(0);  //curr_in_body.get(r);

                            for (int i = 0; i < newInfections; i++) {

                                genotype_curr_per_patch.get(patch_index).add(genotype);

                            }
                            patch_history.logBirth(t_curr);
                            patch_history.logDeath(Double.POSITIVE_INFINITY);
                            patch_history.logGenotype(genotype);
                            patch_history.logParent(colonizing_genotype);
                            patch_history.logPatch(patch_index);

                            List<Integer> new_prevalence = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / tau)]));
                            Collections.fill(new_prevalence, 0);
                            new_prevalence.set(t, newInfections);
                            patch_history.prevalence.add(new_prevalence);
                            genotype++;

                            emptyPatchList.remove((Integer) patch_index);
                            occupyPatchList.add(patch_index);

                            timeRefactory_curr_i[patch_index] = 0.0;

                            j++;
                        }
                        break;
                }
            }

            betweenPatchRates.clear();

            for (int p = 0; p < params.Npatches; p++) {

                List<Integer> genotype_curr = new ArrayList<>();
                genotype_curr.addAll(genotype_curr_per_patch.get(p));

                int Y_temp = genotype_curr.size();

                double infection_rate = r_i*Y_temp;

                if(Y_temp == 0) continue;

                poisson = new Poisson(infection_rate*tau, params.randomGenerator);
                int num = poisson.nextInt();


                int min = Math.min((N - Y_temp), num);


                Integer chosen_genotype = genotype_curr.get(0);
                int genotype_index = patch_history.genotype.indexOf(chosen_genotype);

                int prevalence_curr = patch_history.prevalence.get(genotype_index).get(t - 1);
                int prevalence_next = patch_history.prevalence.get(genotype_index).get(t);


                if (prevalence_next > 0) {

                    prevalence_curr = prevalence_next;
                }

                List<Integer> temp_list = new ArrayList<>(Arrays.asList(new Integer[min]));
                Collections.fill(temp_list, chosen_genotype);
                genotype_curr.addAll(temp_list);
                patch_history.prevalence.get(genotype_index).set(t, prevalence_curr + min);

                if (genotype_curr.size() >= threshold && !thresholdList.contains(new Integer(p))) {

                    thresholdList.add(p);
                }

                genotype_curr_per_patch.replace(p, genotype_curr);
                withinPatchRates.clear();
            }

            curr_in_body.clear();

        }

        if(writeOutput) {
            try {
                writer1.close();
                writer2.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }


    private void initialiseHistory(infectionHistory history, int Y_init, double t_max) {

        int i = 0;
//        while(i < n_patches) {

//            if(Y_curr_i[i] > 0) {
        history.logGenotype(genotype);
        history.logBirth(Double.NEGATIVE_INFINITY);
        history.logDeath(Double.POSITIVE_INFINITY);
        history.logParent((int) Double.NEGATIVE_INFINITY);
        List<Integer> prevalence_through_time = new ArrayList<>(Arrays.asList(new Integer[(int) (t_max / params.tau)]));
        Collections.fill(prevalence_through_time, 0);
        prevalence_through_time.set(0, Y_init);
        history.prevalence.add(prevalence_through_time);
        history.logPatch(i);
        genotype++;
//            }
//            i++;
//        }
    }

    public infectionHistory getInfectionHistory() {

        return this.patch_history;
    }


    private void writeOutput(FileWriter writer, double t, int total_infected, int genotypes, int n_patches, List<Double> diversity) {

        try {
            writer.write(t+"\t"+total_infected+"\t"+genotypes+"\t"+n_patches+"\t"+diversity.get(0)+"\t"+diversity.get(1) +"\n");
            writer.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private int commonAncestor(Integer i1, Integer i2) {

        Integer commonAnc = (int)Double.NEGATIVE_INFINITY;
        Integer lineageA = i1;
        Integer lineageB = i2;


        Set<Integer> ancestry = new HashSet<Integer>();
        while (true) {

            Integer lineage_a = patch_history.getParent(lineageA);
            Integer lineage_b = patch_history.getParent(lineageB);


            if (lineage_a.intValue() != (int) Double.NEGATIVE_INFINITY) {

                lineageA = lineage_a;
                if (!ancestry.add(lineageA)) {
                    commonAnc = lineageA;
                    break;
                }
            }
            if (lineage_b.intValue() != (int) Double.NEGATIVE_INFINITY) {

                lineageB = lineage_b;
                if (!ancestry.add(lineageB)) {
                    commonAnc = lineageB;
                    break;
                }
            }
            if (lineage_a.intValue() == (int) Double.NEGATIVE_INFINITY && lineage_b.intValue() == (int) Double.NEGATIVE_INFINITY) {
                break;
            }

        }


        return commonAnc;
    }

    private double distance(Integer i1, Integer i2) {

        Integer ancestor = commonAncestor(i1,i2);
        if(ancestor.intValue() != (int)Double.NEGATIVE_INFINITY) {

            double ancestorDistance= patch_history.getBirth(ancestor);

            if(Double.isInfinite(ancestorDistance)) {
                ancestorDistance = 0;
            }

            double distA = patch_history.getBirth(i1) - ancestorDistance;
            double distB = patch_history.getBirth(i2) - ancestorDistance;
            return distA + distB;
        }
        else{
            return 0;

        }
    }

    private double geneticDistance(Integer i1, Integer i2) {

        Integer ancestor = commonAncestor(i1,i2);

        double distA = patch_history.getMutationsFromParent(i1, ancestor);
        double distB = patch_history.getMutationsFromParent(i2, ancestor);
        return distA + distB;

    }

    public void updatePatchStatistics(FileWriter writer, double t_curr, boolean writeOutput) {

//        total_infected = 0;
//        total_genotypes = 0;

        for (Map.Entry<Integer, List<Integer>> entry : genotype_curr_per_patch.entrySet()) {

//            for (Integer patch: genotype_curr_per_patch.keySet()) {
            Integer patch = entry.getKey();
            List<Integer> genotype_curr = entry.getValue();

            int Y_i = genotype_curr.size();
            int unique_genotypes = new HashSet<>(genotype_curr).size();
//            total_infected += Y_i;
//            total_genotypes += unique_genotypes;
            int write = -1;

            if (writeOutput == false) {
                continue;
            } else {
                write = 0;
            }

            int check = Y_i + write;


            if (check > 0) {

                try {
                    writer.write(t_curr + "\t" + (patch + 1) + "\t" + Y_i +
                            "\t" + unique_genotypes + "\n");
                    writer.flush();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }

    }

    public void update_curr_in_body(List<Integer> curr_in_body, int t) {
        for (int i = 0; i < patch_history.prevalence.size(); i++) {

            //List<Integer> prevalence = patch_history.prevalence.get(i);
            int prev_curr = patch_history.prevalence.get(i).get(t - 1);
            patch_history.prevalence.get(i).set(t, prev_curr);


            Integer genotype_i = patch_history.genotype.get(i);

            if (prev_curr > 0) {

                List<Integer> tmp = Collections.nCopies(prev_curr, genotype_i);
                curr_in_body.addAll(tmp);
            }

        }

    }

    public List<Double> updateDiversity(boolean getDivergence) {

        double diversity = 0.0;
        double tmrca1 = 0.0;
        int sampleCount1 = 0;

        //double geneticDiversity = 0.0;
        double divergence = 0.0;



            List<Integer> weights = new ArrayList<>();
            for(Integer p: occupyPatchList) {

                weights.add(genotype_curr_per_patch.get(p).size());

            }

            List<Integer> chosenPatches = utils.multinomSamp(weights, occupyPatchList, 100);


            List<Integer> genotype_curr = new ArrayList<>();
            chosenPatches.forEach(i -> genotype_curr.add(genotype_curr_per_patch.get(i).get(0)));






        List<Double> results = new ArrayList<>();

        if(genotype_curr.size()==0) {
            results.add(diversity);
            results.add(tmrca1);
            //results.add(geneticDiversity);
            if(getDivergence) {
                results.add(divergence);
            }

            return results;
        }


        for (int i = 0; i < genotype_curr.size(); i++) {

            Integer vA = genotype_curr.get(Uniform.staticNextIntFromTo(0, genotype_curr.size()-1));
            Integer vB = genotype_curr.get(Uniform.staticNextIntFromTo(0, genotype_curr.size()-1));

            //if (vA.intValue() != (int)Double.NEGATIVE_INFINITY && vB.intValue() != (int)Double.NEGATIVE_INFINITY) {
            double dist = distance(vA, vB);
            //double geneticDist = geneticDistance(vA, vB);
            double div_a = patch_history.getBirth(vA);
            double div_b = patch_history.getBirth(vB);

            if(Double.isInfinite(div_a)) {
                div_a = 0.0;
            }
            if(Double.isInfinite(div_b)) {
                div_b = 0.0;
            }

            diversity += dist;
            //geneticDiversity += geneticDist;
            divergence += div_a;
            divergence += div_b;


            if (dist > tmrca1) {
                tmrca1 = dist;
            }
            sampleCount1 += 1;
        }

        if (sampleCount1 > 0) {
            diversity /= sampleCount1;
            //geneticDiversity /= (double) sampleCount1;
            divergence /= (double) sampleCount1*2;



        }

        tmrca1 /= 2.0;

        results.add(diversity);
        results.add(tmrca1);
        //results.add(geneticDiversity);
        if(getDivergence) {
            results.add(divergence);
        }


        return results;

    }

    public int choosePatch(double[] prob) {

        double x = Uniform.staticNextDouble();

        int index = 0;
        double i = prob[index];
        int patch = occupyPatchList.get(index);


        while(i < x || genotype_curr_per_patch.get(patch).size() == 0) {
            //System.out.println(parent_fitness+","+x);
            index++;
            if(index == (prob.length-1)) {
                break;
            }
            i = prob[index];
            patch = occupyPatchList.get(index);
            //x = Uniform.staticNextDouble();
        }

        //Double real_fitness = iMatrix.getFitness(parent_index);
        //System.out.println(patch+", "+genotype_curr_per_patch.get(patch).size());
        return occupyPatchList.get(index);
    }

    public static void main(String [] args) {

        patchSim test = new patchSim();
        params inputParams = new params();
        inputParams.runTime = Double.parseDouble(args[0]);
        inputParams.startTime = Double.parseDouble(args[1]);
        inputParams.Npatches = Integer.parseInt(args[2]);
        inputParams.S = Integer.parseInt(args[3]);
        inputParams.I = Integer.parseInt(args[4]);
        inputParams.nu = Double.parseDouble(args[5]);
        inputParams.c = Double.parseDouble(args[6]);
        inputParams.r = Double.parseDouble(args[7]);

        test.runSim(inputParams);

    }


    private  Comparator<Integer>  comparator1 = new Comparator<Integer>() {
        @Override
        public int compare(Integer o1, Integer o2) {


            double diff = genotype_curr_per_patch.get(o1).size() - genotype_curr_per_patch.get(o2).size();
            if(diff == 0) {

                return 0;
            }
            else if(diff > 0) {
                return -1;
            }
            else{
                return 1;
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



}
