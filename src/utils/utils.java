package utils;

import cern.jet.random.Uniform;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class utils {


    public static List<Integer> multinomSamp(List<Integer> weights, List<Integer> lineages, int N) {


        List<Integer> sample = new ArrayList<>();

        int sum = weights.stream().mapToInt(Integer::intValue).sum();

        double[] pdf = new double[weights.size()];

        List<Integer> weights_copy = new ArrayList<>();
        weights_copy.addAll(weights);

        weights_copy.forEach(s -> {
            pdf[weights_copy.indexOf(s)] = s / (double) sum;
            weights_copy.set(weights_copy.indexOf(s), -1);
        });


        for (int i = 1; i < pdf.length; i++) {

            pdf[i] = pdf[i] + pdf[i - 1];

        }

        for (int i = 0; i < N; i++) {

            double random = Uniform.staticNextDouble();

            int count = 0;
            Integer chosen = lineages.get(count);
            while (random >= pdf[count]) {

                count++;
                if (count == lineages.size()) {
                    chosen = lineages.get(lineages.size() - 1);
                    break;
                }
                chosen = lineages.get(count);

            }
            sample.add(chosen);

        }

        return sample;
    }

    public static Integer sample(List<Integer> weights, List<Integer> list) {

        double sum = weights.stream().mapToInt(Integer::intValue).sum();

        double[] pdf = new double[weights.size()];

        List<Integer> weights_copy = new ArrayList<>();
        weights_copy.addAll(weights);

        weights_copy.forEach(s -> {
            pdf[weights_copy.indexOf(s)] = s / (double) sum;
            weights_copy.set(weights_copy.indexOf(s), -1);
        });


        for (int i = 1; i < pdf.length; i++) {

            pdf[i] = pdf[i] + pdf[i - 1];

        }

        double random = Uniform.staticNextDouble();

        int count = 0;
        Integer chosen = list.get(count);
        while (random >= pdf[count]) {

            count++;
            if (count == list.size()) {
                chosen = list.get(list.size() - 1);
                break;
            }
            chosen = list.get(count);

        }

        return chosen;
    }

    public static int sample(List<Double> weights) {

        double sum = weights.stream().mapToDouble(Double::doubleValue).sum();

        double[] pdf = new double[weights.size()];

        List<Double> weights_copy = new ArrayList<>();
        weights_copy.addAll(weights);

        weights_copy.forEach(s -> {
            pdf[weights_copy.indexOf(s)] = s / (double) sum;
            weights_copy.set(weights_copy.indexOf(s), -1.);
        });


        for (int i = 1; i < pdf.length; i++) {

            pdf[i] = pdf[i] + pdf[i - 1];

        }

        double random = Uniform.staticNextDouble();

        int count = 0;
        //int chosen = list.get(count);
        while (random >= pdf[count]) {

            count++;
            if (count == weights.size()) {
                count = weights.size() - 1;
                break;
            }

        }

        return count;
    }

    public static int sampleIntegerWeights(List<Integer> weights) {

        double sum = weights.stream().mapToDouble(Integer::intValue).sum();

        double[] pdf = new double[weights.size()];

        List<Integer> weights_copy = new ArrayList<>();
        weights_copy.addAll(weights);

        weights_copy.forEach(s -> {
            pdf[weights_copy.indexOf(s)] = s / (double) sum;
            weights_copy.set(weights_copy.indexOf(s), -1);
        });


        for (int i = 1; i < pdf.length; i++) {

            pdf[i] = pdf[i] + pdf[i - 1];

        }

        double random = Uniform.staticNextDouble();

        int count = 0;
        //int chosen = list.get(count);
        while (random >= pdf[count]) {

            count++;
            if (count == weights.size()) {
                count = weights.size() - 1;
                break;
            }

        }

        return count;
    }



    public static List<Integer> intersect(List<Integer> list1, List<Integer> list2) {


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

    public static List<Integer> find(List<Integer> list, int value) {


        List<Integer> indices = IntStream.range(0, list.size())
                .filter(i -> list.get(i) == value)
                .boxed()
                .collect(Collectors.toList());


        return indices;

    }

    public List<Integer> getIndices(List list) {

        List<Integer> indices = IntStream.range(0, list.size())
                .boxed()
                .collect(Collectors.toList());

        return indices;
    }

    public int sumArray(int[] array) {

        int i = 0;
        int sum = 0;
        while (i < array.length) {

            sum += array[i];
            i++;
        }
        return sum;
    }


    public int chooseBin(double[] prob, double randomNo) {

        int bin = 0;

        while (prob[bin] < randomNo) {

            bin++;

            if (bin == prob.length) {
                break;
            }

        }

        return bin;

    }

    private List<Integer> deepCopy(List<Integer> list) {
        List<Integer> copy = new ArrayList<Integer>(list.size());
        for (Integer element : list) {
            Integer elementCopy = element;
            copy.add(elementCopy);
        }
        return copy;
    }


}