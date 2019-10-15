import cern.jet.random.Binomial;
import cern.jet.random.Empirical;
import cern.jet.random.EmpiricalWalker;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.cert.CollectionCertStoreParameters;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.LongStream;


public class extractGenealogy {

    List<String> names;
    sampledLineages lineages;
    Map<Integer, List<Integer>> parentLineages;
    int[][] b;
    double[] d;
    Map<Integer, String> nodeNames;
    Map<Integer, Double> nodeHeights;

    public extractGenealogy() {

        names = new ArrayList<>();
        nodeNames = new HashMap<>();
        nodeHeights = new HashMap<>();

    }

    public void sampleLineages() {

        //sample lineages from the simulation from which we will construct a transmission tree

    }

    public void getGenealogy(infectionHistory history, Map<Integer, List<Integer>> N_per_patch_over_time, double tau) {

        List<Integer> sampledLineages = new ArrayList<>();
        List<Double> sampleTimes = new ArrayList<>();

        sampledLineages.addAll(lineages.getSampledLineages());
        sampleTimes.addAll(lineages.getSampledTimes());

//        Map<Integer, String> nodeNames = new HashMap<>();
//        Map<Integer, Double> nodeHeights = new HashMap<>();
        Map<Integer, Integer> genotypeMap = new HashMap<>();

        List<String> names = createNames(history);

        //int n_lineages = sampledLineages.size();

        fillNodeMaps(nodeNames, nodeHeights, genotypeMap);

        List<Integer> genotypes = new ArrayList<>(new HashSet<>(sampledLineages));
        List<Double> genotypeAge = new ArrayList<>();


        Collections.sort(genotypes);
        genotypes.forEach(i -> genotypeAge.add(history.getBirth(history.genotype.indexOf(i))));


        genotypes.sort(Comparator.comparing((Integer i) -> genotypeAge.get(genotypes.indexOf(i))));
        //genotypes.sort(Comparator.comparingDouble(genotypes::indexOf));

                //(int)(genotypeAge.get(genotypes.indexOf(i2))-genotypeAge.get(genotypes.indexOf(i1))));

        genotypes.forEach(g -> System.out.println(g+ " : "+history.getBirth(g)));


//        Map<Integer, Integer> genotypeFrequency = new HashMap<>();
//        genotypes.forEach(g -> genotypeFrequency.put(g, Collections.frequency(sampledLineages, g)));


        genotypes.forEach(g -> {

            List<Integer> nSampled = new ArrayList<>();
            List<Double> times = new ArrayList<>();

            IntStream.range(0, sampledLineages.size()).
                    forEach(i -> {
                        if(g.equals(sampledLineages.get(i))) {
                            times.add(sampleTimes.get(i));
                        }
                    });


            //nSampled_per_genotype.put(g, times);

            //double maxTime = Collections.max(times);
            //IntStream.range(0, times.size()).forEach(i -> times.set(i, maxTime-times.get(i)));
            List<Double> sample_times = new ArrayList<>(new HashSet<>(times));
            Collections.sort(sample_times, Comparator.comparingDouble(Double::intValue).reversed());
            //System.out.println(sample_times);


            sample_times.forEach(st -> nSampled.add(times.stream()
                    .filter(j -> j.equals(st))
                    .collect(Collectors.toList())
                    .size()));




            if(nSampled.get(0) > 1 || sample_times.size() > 1) {
                System.out.println(g+" : "+ nSampled);
                System.out.println(sample_times);
                double maxTime = Collections.max(sample_times);
                IntStream.range(0, sample_times.size()).forEach(i -> sample_times.set(i, params.runTime-sample_times.get(i)));
                System.out.println(sample_times);

                int index = history.genotype.indexOf(g);
                int patch = history.getPatch(index);
                coalSim sim = new coalSim(N_per_patch_over_time.get(patch), tau, history);
                //N_per_patch_over_time.get(patch), tau, nodeNames, nodeHeights);
                sim.sampleTree(sample_times, nSampled, nodeNames, nodeHeights, genotypeMap, sampledLineages, sampleTimes, g, patch);
                //genotypeFrequency.replace(g, 1);

            }



        });
        List<Integer> completeIndividuals = new ArrayList<>(sampledLineages);
        List<Double> completeSeqTimes = new ArrayList<>(sampleTimes);

        sampledLineages.removeAll(new ArrayList<>(Arrays.asList(-1)));
        sampleTimes.removeAll(new ArrayList<>(Arrays.asList(-1.0)));



        List<Integer> individualsSubList = new ArrayList<>(sampledLineages);
        List<Double> timesSubList = new ArrayList<>(sampleTimes);


        int n_lineages = sampledLineages.size();


        b = new int[(n_lineages - 1)][2];
        d = new double[2 * names.size() - 1];


        List<Integer> currIndividuals = new ArrayList<>(sampledLineages);
        List<Double> currSeqTimes = new ArrayList<>(sampleTimes);

        // I think this needs to come earlier in the code...2019/10/08


        boolean condition = true;

        int loc_b = 0;
        int x_event = 1;

        List<Integer> coalescedList = new ArrayList<>();

        parentLineages = getParentLineages(currIndividuals, history);

        while (condition) {

            System.out.println(".coal event " + x_event);

            coalescentEvent coalescentEvent = findMostRecentCoalescence(currIndividuals, history);

            double timeOfCoalescence = coalescentEvent.getTimeToCoalescence();
            Integer coal_parent = coalescentEvent.getCoalParent();
            int[] coalDaughters = coalescentEvent.getCoalDaughters();

            if (timeOfCoalescence == Double.NEGATIVE_INFINITY) {
                break;
            }

            List<Integer> indiv1_index = find(completeIndividuals, coalDaughters[0]);
            List<Integer> indiv2_index = find(completeIndividuals, coalDaughters[1]);

            List<Integer> non_coalescedList1 = new ArrayList<>();
            List<Integer> non_coalescedList2 = new ArrayList<>();

            for (Integer i : indiv1_index) {

                Integer x = i + 1;
                if (!coalescedList.contains(x)) {

                    non_coalescedList1.add(i);
                }
            }

            for (Integer i : indiv2_index) {

                Integer x = i + 1;
                if (!coalescedList.contains(x)) {

                    non_coalescedList2.add(i);
                }
            }


            int index1 = non_coalescedList1.get(0);
            int index2 = non_coalescedList2.get(0);

            if (index1 == index2) {

                index2 = non_coalescedList1.get(1);

            }


            int d_index1 = individualsSubList.indexOf(coalDaughters[0]);
            int d_index2 = individualsSubList.indexOf(coalDaughters[1]);

            d[index1] = Math.abs(completeSeqTimes.get(index1) - timeOfCoalescence);
            d[index2] = Math.abs(completeSeqTimes.get(index2) - timeOfCoalescence);


            int[] children = new int[2];

            b[loc_b] = new int[2];
            children[0] = index1 + 1;
            children[1] = index2 + 1;
            if (coalescedList.contains(children[0]) || coalescedList.contains(children[1])) {

                System.out.println("coalesced lineage already");


            }
            coalescedList.add(children[0]);
            coalescedList.add(children[1]);

            Arrays.sort(children);
            b[loc_b] = children;

            int history_index = history.genotype.indexOf(coal_parent);
            String nodeName = "'coal_" + x_event + "_patch_" + (history.getPatch(history_index)+1) +
                    "_node_" + (loc_b + 1) + "_" + coal_parent + "_" + (timeOfCoalescence) + "'";


//            String name_1 = names.get(children[0] - 1);
//            String name_2 = names.get(children[1] - 1);
//
//            names.set(children[0] - 1, name_1); //new_name1);
//            names.set(children[1] - 1, name_2); //new_name2);

            names.add(nodeName);

            loc_b += 1;

            completeIndividuals.set(index1, -1);
            completeIndividuals.set(index2, -1);

            completeIndividuals.add(coal_parent);

            completeSeqTimes.add(timeOfCoalescence);

            currIndividuals.remove((Integer) coalDaughters[0]);
            currIndividuals.remove((Integer) coalDaughters[1]);

            currIndividuals.add(coal_parent);

            individualsSubList.add(coal_parent);
            timesSubList.add(timeOfCoalescence);
            individualsSubList.set(d_index1, -1);
            individualsSubList.set(d_index2, -1);



            n_lineages = currIndividuals.size();

            parentLineages = getParentLineages(currIndividuals, history);

            x_event++;

            if (n_lineages == 1) {
                condition = false;
            }
        }

        System.out.println(n_lineages + " have not coalesced");

        while (condition) {
            if (n_lineages == 1) {
                break;
            }

            List<Integer> indices = IntStream.range(0, currIndividuals.size()).
                    boxed().
                    collect(Collectors.toList());  //getIndices(currIndividuals);
            //Collections.shuffle(indices);

            // these should be unique
            int index_1 = indices.get((int) Uniform.staticNextDouble() * indices.size());
            indices.remove((Integer) index_1);
            int index_2 = indices.get((int) Uniform.staticNextDouble() * indices.size());

            int[] coalDaughters = new int[2];
            //coalDaughters might be the same
            coalDaughters[0] = currIndividuals.get(index_1);
            coalDaughters[1] = currIndividuals.get(index_2);

            List<Integer> indiv1_index = find(completeIndividuals, coalDaughters[0]);
            List<Integer> indiv2_index = find(completeIndividuals, coalDaughters[1]);

            List<Integer> non_coalescedList1 = new ArrayList<Integer>();
            List<Integer> non_coalescedList2 = new ArrayList<Integer>();

            for (Integer i : indiv1_index) {

                Integer x = i + 1;
                if (!coalescedList.contains(x)) {

                    non_coalescedList1.add(i);
                }
            }

            for (Integer i : indiv2_index) {

                Integer x = i + 1;
                if (!coalescedList.contains(x)) {

                    non_coalescedList2.add(i);
                }
            }

            int index1 = non_coalescedList1.get(0);
            int index2 = non_coalescedList2.get(0);

            //this means non_coalesced lists 1 and 2 should be the same, and should contain at most
            if (index1 == index2) {

                System.out.println("nc list1: " + non_coalescedList1);
                System.out.println("nc list2: " + non_coalescedList2);

                index2 = non_coalescedList1.get(1);//tools.lastIndex(non_coalescedList1));

            }

            int[] children = new int[2];
            b[loc_b] = new int[2];
            children[0] = index1 + 1;
            children[1] = index2 + 1;

            if (children[0] == children[1]) {

                System.out.println("Stop");
            }
            Arrays.sort(children);

            b[loc_b] = children;

            coalescedList.add(children[0]);
            coalescedList.add(children[1]);

            double timeOfCoalescence = 0;

            d[index1] = (completeSeqTimes.get(index1) - timeOfCoalescence);
            d[index2] = (completeSeqTimes.get(index2) - timeOfCoalescence);


            int coal_parent = Collections.max(completeIndividuals) + 1;


            String nodeName = "'coal_" + x_event + "_patch_-1" + "_node_" + (loc_b + 1) + "_" + (timeOfCoalescence) + "'";
//"'coal_"+x_event+"_mutations_-1_rate_-1_patch_-1"+"_node_"+(loc_b+1)+"_"+(timeOfCoalescence)+"'";

            x_event++;

            names.add(nodeName);

            loc_b += 1;


            completeIndividuals.add(coal_parent);
            completeSeqTimes.add(timeOfCoalescence);
            currIndividuals.add(coal_parent);

            currIndividuals.remove((Integer) coalDaughters[0]);

            currIndividuals.remove((Integer) coalDaughters[1]);

            n_lineages = currIndividuals.size();

        }
        System.out.println("names b: " + names.size());


        d[d.length - 1] = -50;


        Collections.sort(coalescedList);


    }

    public Map<Integer, List<Integer>> getParentLineages(List<Integer> sampledLineages, infectionHistory history) {

        //get parents of particular set of lineages of interest

        Map<Integer, List<Integer>> parentLineages = new HashMap();

        for (int i = 0; i < sampledLineages.size(); i++) {

            List<Integer> lineages = getLineage(sampledLineages.get(i), history);
            parentLineages.put(sampledLineages.get(i), lineages);

        }

        return parentLineages;

    }

    public List<Integer> getLineage(Integer indiv_lineage, infectionHistory history) {

        List<Integer> parents = new ArrayList<>();
        parents.addAll(history.parent);

        int i = 0;

        List<Integer> this_lineage = new ArrayList<Integer>();

        this_lineage.add(indiv_lineage);

        boolean condition = true;
        while (condition) {


            Integer lineage_parent = parents.get(history.genotype.indexOf(this_lineage.get(i)));
            if (lineage_parent < 0) {

                return this_lineage;
            } else {

                this_lineage.add(lineage_parent);
                i++;
            }
        }

        return this_lineage;

    }

    public coalescentEvent findMostRecentCoalescence(List<Integer> current_indiv, infectionHistory history) {

        double mostRecentTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;
        double coalTime1 = Double.NEGATIVE_INFINITY;
        double coalTime2 = Double.NEGATIVE_INFINITY;
        int n_individuals_sampled = current_indiv.size();
        int[] coalDaughters = new int[2];
        Integer coalParent = (int) Double.NEGATIVE_INFINITY;
        List<Double> births = history.getCompleteBirths();
        List<Double> deaths = history.getCompleteDeaths();

        coalescentEvent coalescentEvent = null;
        List<Integer> currIndividuals = new ArrayList<Integer>(current_indiv);

        for (int i = 0; i < n_individuals_sampled; i++) {


            for (int j = (i + 1); j < n_individuals_sampled; j++) {

                Integer indiv_i = currIndividuals.get(i);
                Integer indiv_j = currIndividuals.get(j);

                List<Integer> parentsInCommon = intersect(parentLineages.get(indiv_i), parentLineages.get(indiv_j));
                //System.out.println(parentsInCommon);

                if (!parentsInCommon.isEmpty()) {

                    Integer this_parent = Collections.max(parentsInCommon);

                    List<Integer> p1 = parentLineages.get(indiv_i);
                    List<Integer> p2 = parentLineages.get(indiv_j);

                    int loc1 = p1.indexOf(this_parent);
                    int loc2 = p2.indexOf(this_parent);

                    if (currIndividuals.get(i).equals(this_parent)) {

                        int index = history.genotype.indexOf(currIndividuals.get(i));
                        Integer parent_genotype = history.parent.get(index);
                        if(history.genotype.contains(parent_genotype)) {
                            int parent_index = history.genotype.indexOf(parent_genotype);
                            coalTime1 = births.get(parent_index);
                        }
                        else {
                            coalTime1 = deaths.get(index);
                        }


                    } else {

                        int parent_index = history.genotype.indexOf(parentLineages.get(indiv_i).get(loc1 - 1));
                        coalTime1 = births.get(parent_index);//parentLineages.get(indiv_i).get(loc1-1));
                    }

                    if (currIndividuals.get(j).equals(this_parent)) {

                        int index = history.genotype.indexOf(currIndividuals.get(j));
                        Integer parent_genotype = history.parent.get(index);
                        if(history.genotype.contains(parent_genotype)) {
                            int parent_index = history.genotype.indexOf(parent_genotype);
                            coalTime2 = births.get(parent_index);
                        }
                        else {
                            coalTime2 = deaths.get(index);
                        }

                    } else {

                        int parent_index = history.genotype.indexOf(parentLineages.get(indiv_j).get(loc2 - 1));
                        coalTime2 = births.get(parent_index); //parentLineages.get(indiv_j).get(loc2-1));
                    }

                    thisTimeOfCoalescence = Math.min(coalTime1, coalTime2);


                    if (thisTimeOfCoalescence > mostRecentTimeOfCoalescence) {

                        if (!currIndividuals.get(i).equals(currIndividuals.get(j))) {

                            //thisTimeOfCoalescence = Double.NEGATIVE_INFINITY;

                            coalDaughters = new int[2];
                            coalDaughters[0] = currIndividuals.get(i);
                            coalDaughters[1] = currIndividuals.get(j);
                            coalParent = this_parent;
                            mostRecentTimeOfCoalescence = thisTimeOfCoalescence;
                        }
                    }
                }
            }
        }

        coalescentEvent = new coalescentEvent(mostRecentTimeOfCoalescence, coalParent, coalDaughters);

        return coalescentEvent;

    }

    private List<Integer> intersect(List<Integer> list1, List<Integer> list2) {


        Set<Integer> set1 = new HashSet<Integer>();
        Set<Integer> set2 = new HashSet<Integer>();

        set1.addAll(list1);
        set2.addAll(list2);
        //Collections.copy(newList, list1);

        set1.retainAll(set2);

        List<Integer> intersectList = new ArrayList<Integer>();
        intersectList.addAll(set1);
        return intersectList;
    }

    public sampledLineages getSampledGenotypes(infectionHistory history, double startTime) {

        Map<Integer, Double> lineagesMap = new HashMap<>();

        List<Integer> sampledLineages = new ArrayList<>();
        List<Double> sampledTimes = new ArrayList<>();
        Set<Integer> lineagesSet = new HashSet<>();

        int interval = (int) Math.ceil((params.runTime - startTime) / 10);

        double t = startTime;
        while (t < params.runTime) {

            t = t + interval;

            System.out.println(t);

            List<Integer> lineages = new ArrayList<>();

            for (int i = 0; i < history.getCompleteBirths().size(); i++) {

                if (history.birth.get(i) <= t && history.death.get(i) >= t && history.birth.get(i) > startTime) {

                    lineages.add(history.genotype.get(i));
                    lineagesMap.put(history.genotype.get(i), history.death.get(i));
                }
            }

            Collections.shuffle(lineages);
            if (lineages.size() > params.n_samples_per_time) {
                lineagesSet.addAll(lineages.subList(0, params.n_samples_per_time));
            } else {
                lineagesSet.addAll(lineages);
            }

        }

        sampledLineages.addAll(lineagesSet);
        Collections.sort(sampledLineages);
        for (int i = 0; i < sampledLineages.size(); i++) {


            double sampleTime = lineagesMap.get(sampledLineages.get(i));
            if (Double.isInfinite(sampleTime)) {
                sampleTime = params.runTime;
            }
            System.out.println(sampledLineages.get(i) + " " + sampleTime);
            sampledTimes.add(sampleTime);

        }

        return new sampledLineages(sampledLineages, sampledTimes);

    }

    public sampledLineages getTestSampledLineages(infectionHistory history, double startTime, double tau) {

        Map<String, Double> lineagesMap = new HashMap<>();

        List<Integer> sampledLineages = new ArrayList<>();
        List<Double> sampledTimes = new ArrayList<>();

        int interval = (int) Math.ceil((params.runTime - startTime) / (double)params.interval);

        double t = 300;

        List<Integer> n = new ArrayList<>(Arrays.asList(6,1));
        int count = 0;
        while (t < params.runTime) {


            System.out.println(t);

            List<Integer> n_Sampled = new ArrayList<>();
            List<Integer> lineages = new ArrayList<>();

            for (int i = 0; i < history.getCompleteBirths().size(); i++) {

                if (history.getBirth(i) <= t && history.getDeath(i) > t && history.getBirth(i) > startTime) {

                    Integer genotype = history.genotype.get(i);
                    lineages.add(genotype);

                    //lineagesMap.put(genotype+"_"+t, t);
                    n_Sampled.add(history.prevalence.get(i).get((int) (t / tau) - 1));
                    //n_Sampled.add(history.getPrevalenceAtTimet(genotype, (int)(t/tau)));

                    if (history.prevalence.get(i).get((int) (t / tau) - 1) == 0) {
                        System.out.println("> " + genotype + ", " + history.getBirth(i) + ", " + history.getDeath(i) + ", " + t);
                    }

                }
            }

            //System.out.println(params.n_samples_per_time);
            //List<Integer> lineages_t = multinomSamp(n_Sampled, lineages, n.get(count));//
            List<Integer> lineages_t =  new ArrayList<>(Collections.nCopies(2, lineages.get((int)Math.floor(Math.random()*lineages.size()))));
            Collections.sort(lineages_t);
            sampledLineages.addAll(lineages_t);
            sampledTimes.addAll(Collections.nCopies(lineages_t.size(), t));
            count++;
            t = t + 50;

        }


        return new sampledLineages(sampledLineages, sampledTimes);

    }

    // this method samples lineages proportional to their prevalence
    public sampledLineages getSampledLineages(infectionHistory history, double startTime, double tau) {

        Map<String, Double> lineagesMap = new HashMap<>();

        List<Integer> sampledLineages = new ArrayList<>();
        List<Double> sampledTimes = new ArrayList<>();

        int interval = (int) Math.ceil((params.runTime - 200.0) / (double)params.interval);

        double t = 200.0;
        while (t < params.runTime) {

            System.out.println(t);

            List<Integer> n_Sampled = new ArrayList<>();
            List<Integer> lineages = new ArrayList<>();

            for (int i = 0; i < history.getCompleteBirths().size(); i++) {

                if (history.getBirth(i) <= t && history.getDeath(i) > t && history.getBirth(i) > startTime) {

                    Integer genotype = history.genotype.get(i);
                    lineages.add(genotype);

                    //lineagesMap.put(genotype+"_"+t, t);
                    n_Sampled.add(history.prevalence.get(i).get((int) (t / tau) - 1));
                    //n_Sampled.add(history.getPrevalenceAtTimet(genotype, (int)(t/tau)));

                    if (history.prevalence.get(i).get((int) (t / tau) - 1) == 0) {
                        System.out.println("> " + genotype + ", " + history.getBirth(i) + ", " + history.getDeath(i) + ", " + t);
                    }

                }
            }

            //System.out.println(params.n_samples_per_time);
            List<Integer> lineages_t = multinomSamp(n_Sampled, lineages, params.n_samples_per_time);
            Collections.sort(lineages_t);
            sampledLineages.addAll(lineages_t);
            sampledTimes.addAll(Collections.nCopies(lineages_t.size(), t));

            t = t + interval;
        }

//        Collections.sort(sampledLineages);
//        for (int i = 0; i < sampledLineages.size(); i++) {
//
//
//            double sampleTime = lineagesMap.get(sampledLineages.get(i));
//            if (Double.isInfinite(sampleTime)) {
//                sampleTime = params.runTime;
//            }
//            System.out.println(sampledLineages.get(i) + " " + sampleTime);
//            sampledTimes.add(sampleTime);
//
//        }

        return new sampledLineages(sampledLineages, sampledTimes);

    }

    // this method samples unique genotypes is now deprecated
    public sampledLineages getSampledLineages(infectionHistory history, int nLineages) {


        //List<Integer> lineages = new ArrayList<>();
        List<Integer> lineagesList = new ArrayList<>();

        List<Integer> sampledLineages = new ArrayList<>();
        List<Double> sampledTimes = new ArrayList<>();


        int count = 1;
        for (int i = 0; i < history.getCompleteBirths().size(); i++) {

            if (history.death.get(i) >= params.startTime && history.birth.get(i) > 0) {

                lineagesList.add(history.genotype.get(i));

            }
        }

        Collections.shuffle(lineagesList);
        sampledLineages.addAll(lineagesList.subList(0, nLineages));

        for (int i = 0; i < sampledLineages.size(); i++) {

            int index = history.genotype.indexOf(sampledLineages.get(i));

            count++;

            double sampleTime = history.death.get(index);
            if (Double.isInfinite(sampleTime)) {
                sampleTime = params.runTime;
            }
            sampledTimes.add(sampleTime);

        }


        return new sampledLineages(sampledLineages, sampledTimes);

    }

    private void fillNodeMaps(Map<Integer, String> nodeNames, Map<Integer, Double> nodeHeights, Map<Integer, Integer> genotypeMap) {

        List<Integer> sampledLineages = new ArrayList<>();
        List<Double> sampleTimes = new ArrayList<>();
        sampledLineages.addAll(lineages.getSampledLineages());
        sampleTimes.addAll(lineages.getSampledTimes());

        for (int i = 0; i < sampledLineages.size(); i++) {
            nodeNames.put(i + 1, "'" + names.get(i) + "'");
            nodeHeights.put(i + 1, sampleTimes.get(i));
            genotypeMap.put(i + 1, sampledLineages.get(i));
        }
    }

    private List<String> createNames(infectionHistory history) {

        List<Integer> sampledLineages = lineages.getSampledLineages();
        List<Double> sampleTimes = lineages.getSampledTimes();
        for (int i = 0; i < sampledLineages.size(); i++) {

            int index = history.genotype.indexOf(sampledLineages.get(i));
            String name = "sample_" + (i+1) + "_genotype_"+ history.getId(index)+"_patch_" + (history.getPatch(index)+1) + "_" + sampleTimes.get(i);
            names.add(name);
        }

        return names;
    }

    private List<Integer> find(List<Integer> list, int value) {

//        List<Integer> indices = new ArrayList<Integer>();
//
//        List<Integer> copyList = new ArrayList<Integer>();

        List<Integer> indices = IntStream.range(0, list.size())
                .filter(i -> list.get(i)==value)
                .boxed()
                .collect(Collectors.toList());

//        while (copyList.contains(value)) {
//
//            int index = copyList.indexOf(value);
//            indices.add(index);
//            copyList.set(index, (int) Double.NEGATIVE_INFINITY);
//            //copyList.remove(index);
//
//        }

        return indices;

    }

    public List<Integer> getIndices(List list) {

        List<Integer> indices = IntStream.range(0, list.size())
                .boxed()
                .collect(Collectors.toList());

        return indices;
    }

    public void writeNexusTree(int[][] coalescentEvents, double[] nodeTimes, List<String> names, File treeFile) {

//        Map<Integer, String> nodeNames = new HashMap<Integer, String>();
//        Map<Integer, Double> nodeHeights = new HashMap<Integer, Double>();
//        Map<Integer, Integer> genotypeMap = new HashMap<Integer, Integer>();

        int node_key = Collections.max(nodeNames.keySet());

        int n_lineages = lineages.getSampledLineages().size();

        int n_nodes = nodeNames.size();
        //fillNodeMaps(nodeNames, nodeHeights, genotypeMap);

        FileWriter writer4 = null;

        try {
            writer4 = new FileWriter(treeFile);


            writer4.write("#NEXUS\n");
            writer4.write("begin taxa;\n");
            writer4.write("\tdimensions ntax=" + n_lineages + ";\n");
            writer4.write("\ttaxlabels\n");

            int n = 0;
            while (n < n_lineages) {
                writer4.write("\t'" + names.get(n) + "'\n");
                n++;
            }
            writer4.write(";\n");
            writer4.write("end;\n");
            writer4.write("begin trees;\n");
            int k = 0;
            String treeString = "";

            while (k < coalescentEvents.length) {

                int child1 = coalescentEvents[k][0];
                int child2 = coalescentEvents[k][1];


                String node1 = "";
                String node2 = "";

                //nodeNames should contain the children names! Just a check that everything is working as it should be!
                if (nodeNames.containsKey(child1)) {

                    node1 = nodeNames.get(child1);
                }
                if (nodeNames.containsKey(child2)) {

                    node2 = nodeNames.get(child2);
                }



                //System.out.println(partsNode1[1] + "\t" + partsNode2[1] + "\t" + names.get(n_lineages + k));

                String node1_string = "";
                String node2_string = "";

                if (child1 <= n_lineages) {

                    node1_string = "(" + node1 + ":" + nodeTimes[child1 - 1] + "[&parent=" + names.get(n_lineages + k) + "]";// ",patch=\"" + patchNode1 + "\"]";//\",mutations="+mutations1+",rate="+rate1+"]";
                   // node1_string = "(" + node1 + ":" + nodeTimes[child1 - 1];
                } else {

                    node1_string = "(" + node1 + "[&parent=" + names.get(n_lineages + k) + "]";
                    //node1_string = "(" + node1  ;//":"+nodeTimes[child1-1]+")"; //+"[&parent="+names.get(n_lineages+k)+",patch=\""+patchNode1+"\",mutations="+mutations1+",rate="+rate1+"]";

                }
                if (child2 <= n_lineages) {

                    node2_string = node2 + ":" + nodeTimes[child2 - 1] + "[&parent=" + names.get(n_lineages + k) + "])"; //+ ",patch=\"" + patchNode2 + "\"])";//\",mutations="+mutations2+",rate="+rate2+"])";
                    //node2_string = node2 + ":" + nodeTimes[child2 - 1] + ")";
                    //")[&parent="+names.get(n_lineages+k)+",patch=\""+patchNode+"\",mutations="+mutations+",rate="+rate+"]";

                } else {
                    node2_string = node2 +"[&parent=" + names.get(n_lineages + k) + "])";//":"+nodeTimes[child2-1]+")" ;//+"[&parent="+names.get(n_lineages+k)+",patch=\""+patchNode2+"\",mutations="+mutations2+",rate="+rate2+"]";

                }

//                treeString = "("+node1+"[&parent="+names.get(n_lineages+k)+",patch=\""+patchNode1+"\",mutations="+mutations1+",rate="+rate1+"]:"+nodeTimes[child1-1]+
//                        ","+node2+"[&parent="+names.get(n_lineages+k)+",patch=\""+patchNode2+"\",mutations="+mutations2+",rate="+rate2+"]:"+nodeTimes[child2-1]
//                        +")[&parent="+names.get(n_lineages+k)+",patch=\""+patchNode+"\",mutations="+mutations+",rate="+rate+"]";

                treeString = node1_string + "," + node2_string + ":" + nodeTimes[n_lineages + k];//"[&parent=" + names.get(n_lineages + k) +"]";//+ ",patch=\"" + patchNode + "\"]";//\",mutations="+mutations+",rate="+rate+"]";
                ;
                //System.out.println(nodeString);


                //treeString = "("+node1+":"+nodeTimes[child1-1]+","+node2+":"+nodeTimes[child2-1]+")";

                System.out.println(treeString);
                nodeNames.remove(child1);
                nodeNames.remove(child2);

                nodeNames.put(node_key + k + 1, treeString);
                nodeHeights.put(node_key + k + 1, nodeTimes[n_lineages + k - 1]);

                k++;
            }

            //add the root node info:
            //treeString = treeString+"[&parent='root'"

            System.out.println("******");
            System.out.println(treeString);
            writer4.write("\t tree TREE1 = [&R] [&R=true]" + treeString + ";\n");
            writer4.write("end;\n");
            writer4.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    class coalescentEvent {

        double timeToCoalescence;
        int coalParent;
        int[] coalDaughters;

        coalescentEvent(double timeCoalescence, Integer coalParent, int[] coalDaughters) {

            this.timeToCoalescence = timeCoalescence;
            this.coalParent = coalParent;
            this.coalDaughters = coalDaughters;
        }

        private double getTimeToCoalescence() {

            return this.timeToCoalescence;
        }

        private Integer getCoalParent() {

            return this.coalParent;
        }

        private int[] getCoalDaughters() {

            return this.coalDaughters;
        }
    }

    class sampledLineages {

        List<Integer> lineages;
        List<Double> times;

        sampledLineages(List<Integer> lineages, List<Double> times) {

            this.lineages = lineages;
            this.times = times;

        }

        private List<Integer> getSampledLineages() {

            return this.lineages;
        }

        private List<Double> getSampledTimes() {

            return this.times;
        }

    }

    public List<Integer> multinomSamp(List<Integer> weights, List<Integer> lineages, int N) {

//        System.out.println(weights);
//        System.out.println(lineages);

        List<Integer> sample = new ArrayList<>();

        int sum = weights.stream().mapToInt(Integer::intValue).sum();

        double[] pdf = new double[weights.size()];

        weights.forEach(s -> {
            pdf[weights.indexOf(s)] = s / (double) sum;
            weights.set(weights.indexOf(s), -1);
        });


        for(int i=1; i < pdf.length; i++) {

            pdf[i] = pdf[i]+pdf[i-1];

        }

        for(int i=0; i < N; i++) {

            double random = Uniform.staticNextDouble();

            int count = 0;
            Integer chosen = lineages.get(count);
            while(random >= pdf[count]) {

                count++;
                chosen = lineages.get(count);

            }
            sample.add(chosen);

        }

        return sample;
    }


    public static void main(String[] args) {



        params inputParams = new params();
        inputParams.runTime = Double.parseDouble(args[0]);
        inputParams.startTime = Double.parseDouble(args[1]);

        inputParams.Npatches = Integer.parseInt(args[2]);
        inputParams.S = Integer.parseInt(args[3]);

        inputParams.I = Integer.parseInt(args[4]);
        inputParams.nu = Double.parseDouble(args[5]); //extinction

        inputParams.c = Double.parseDouble(args[6]); //colonization
        inputParams.U = Double.parseDouble(args[7]);

        inputParams.n_samples_per_time = (int)Math.ceil(50/(double)inputParams.interval);
//
//
        patchSim test = new patchSim();

        test.runSim(inputParams);
        infectionHistory infectionHistory =  test.getInfectionHistory();
        double tau = test.tau;
        extractGenealogy genealogy = new extractGenealogy();

        genealogy.lineages = genealogy.getSampledLineages(infectionHistory, inputParams.startTime, tau);
        //test data
        //genealogy.lineages = genealogy.getTestSampledLineages(infectionHistory, inputParams.startTime, tau);


        //genealogy.lineages =  genealogy.getSampledLineages(infectionHistory, 200);

        List<Integer> sampledLineages = new ArrayList<>();
        sampledLineages.addAll(genealogy.lineages.getSampledLineages());
        System.out.println(sampledLineages);

        //genealogy.parentLineages = genealogy.getParentLineages(sampledLineages, infectionHistory);

        genealogy.getGenealogy(infectionHistory, test.N_per_patch_over_time, tau);
//
        genealogy.writeNexusTree(genealogy.b, genealogy.d, genealogy.names,
                new File("tree_ext_"+ params.nu +"_col_"+ params.c +"_Npatches_"+ params.Npatches +"_patchSize_"+(params.S + params.I)+"_nlineages_"+".tre"));



    }
}
