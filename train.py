#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, random, math, time, sys
import multiprocessing as mul
from itertools import combinations
import cPickle

import const
from util import check_contain_chinese, del_file, copy_file
from libsvm.checkdata import check_data
from libsvm.python.svmutil import *
from libsvm.python.plotroc import *

def data_preprocess(filename, path_name):
    """Do data preprocessing: check if the file is LIBSVM format;
                              tranform the format of file generated from webserver.
    :param filename: file to be parsed.
    """
    line_list = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) == 0 or ((not line[0].isdigit()) and line[0] not in ['+', '-']):
                #print line.strip()[0].isdigit()
                continue
            else:
                line_list.append(line)
    root, extension = os.path.splitext(os.path.basename(filename))
    new_filename = path_name + root + '_new' + extension
    while(os.path.isfile(new_filename)):
        random.seed(time.time())
        new_filename = path_name + root + str(int(random.random() * 100000)) + extension
    with open(new_filename, 'w') as f:
        for i in line_list:
            f.write(i)
            f.write('\n')
    if check_data(new_filename) == 0:
        return new_filename
    elif check_data(new_filename) == 1:
        return False


def trans_labels(file_list):
    """For binary classification, transform the labels to +1 and -1.
       For multiclass classification, transform the labels to 0,1,2,...
    """
    if len(file_list) == 2:
        for index, file in enumerate(file_list):
            new_list = []
            if index == 0:
                with open(file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('+1'):
                            new_list.append(line)
                        else:
                            if line[0].isdigit():
                                new_line = line.replace(line[0], '+1', 1)
                            else:
                                new_line = line.replace(line[:2], '+1', 1)
                            new_list.append(new_line)
                with open(file, 'w') as f:
                    for line in new_list:
                        f.write(line)
                        f.write('\n')
            if index == 1:
                with open(file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('-1'):
                            new_list.append(line)
                        else:
                            if line[0].isdigit():
                                new_line = line.replace(line[0], '-1', 1)
                            else:
                                new_line = line.replace(line[:2], '-1', 1)
                            new_list.append(new_line)
                with open(file, 'w') as f:
                    for line in new_list:
                        f.write(line)
                        f.write('\n')
    elif len(file_list) > 2:
        for index, file in enumerate(file_list):
            new_list = []
            with open(file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(str(index)):
                        new_list.append(line)
                    else:
                        if line[0].isdigit():
                            new_line = line.replace(line[0], str(index), 1)
                        else:
                            new_line = line.replace(line[:2], str(index), 1)
                        new_list.append(new_line)
            with open(file, 'w') as f:
                for line in new_list:
                    f.write(line)
                    f.write('\n')



def param_selection(file_list, metric, svm_params, process_num, c_range, g_range, bi_or_multi):
    """This is the main process of parameter selection.
    :param pos_seq: positive sequence list.
    :param neg_seq: negative sequence list.
    :param alphabet: list of alphabet.
    :param metric: the performance measurement used for parameter selection.
    :param svm_params: the parameters of svm_train in libsvm.
    """
    start_time = time.time()
    mgr = mul.Manager()
    best_rate = mul.Value('d', 0.0)
    best_params = mgr.list()
    t_v_list = mgr.list()  #used to save the train and validation datasets which have the largest rate

    path = const.GEN_FILE_PATH
    current_path = os.path.dirname(os.path.realpath(__file__))

    # for i in c_range:
    #     print i
    # print 'sep'
    # for j in g_range:
    #     print j

    # pos_file_name = 'pos_svm_file.txt'
    # neg_file_name = 'neg_svm_file.txt'
    process_list = []
    lock = mul.Lock()
    #process_num = mul.cpu_count()
    semph = mul.Semaphore(process_num)

    if len(c_range) <= len(g_range):
        #print 'test1'
        for c in c_range:
            for g in g_range:
                write_path = current_path + path + str(c) + '_' + str(g) + '/'
                if not os.path.exists(write_path):
                    os.makedirs(write_path)
                process_list.append(mul.Process(target=params_selection_one, args=(file_list,
                    metric, write_path, svm_params, lock, semph, best_rate, c, g, best_params, t_v_list, bi_or_multi)))
        for process in process_list:
            process.start()
        for process in process_list:
            process.join()
    elif len(c_range) > len(g_range):
        #print 'test2'
        for g in g_range:
            for c in c_range:
                write_path = current_path + path + str(c) + '_' + str(g) + '/'
                if not os.path.exists(write_path):
                    os.makedirs(write_path)
                process_list.append(mul.Process(target=params_selection_one, args=(file_list,
                    metric, write_path, svm_params, lock, semph, best_rate, c, g, best_params, t_v_list, bi_or_multi)))
        for process in process_list:
            process.start()
        for process in process_list:
            process.join()

    time_cost = time.time() - start_time


    file_path = current_path + path + str(best_params[0]) +'_' + str(best_params[1]) + '/'
    file_name = file_path + 'training_validation_set.txt'
    if os.path.exists(file_path):
        with open(file_name, 'w') as f:
            for part in t_v_list:
                for elem in part:
                    f.write(elem)
                    f.write('\n')

    #train the model using the training and validation sets.
    test_file = file_path + 'test_set.txt'
    if os.path.isfile(test_file):
        evals = one_svm_process(file_name, test_file, metric, svm_params, bi_or_multi)

    print 'The time cost for parameter selection is %.2fs' % time_cost
    return best_params, evals, file_path


def params_selection_one(file_list,  metric, write_path, svm_params, lock, semph, best_rate, c, g, best_params, t_v_list, bi_or_multi):
    """Parameter selection process for one group of parameters: c & g.
    :param file_list: input files list.
    :param metric: the performance measurement used for parameter selection.
    :param write_path: The path for saving the generated files.
    :param svm_params: the parameters of svm_train in libsvm.
    :param lock: The lock used for multiprocessing.
    :param semph: The semaphore used for multiprocessing.
    :param best_rate: Used for storing the best metric.
    :param best_params: Used for storing the best parameters.
    :param t_v_list: Used for storing the training and validation datasets
    generated by the best parameters.
    """
    semph.acquire()
    # pos_file_name = 'pos_svm_file.txt'
    # neg_file_name = 'neg_svm_file.txt'
    # write_file_pos = write_path + pos_file_name
    # write_file_neg = write_path + neg_file_name
    #pos_dataset = feature_ext(pos_seq, '+1', None, alphabet, write_file_pos, k, w, lamada)
    #neg_dataset = feature_ext(neg_seq, '-1', None, alphabet, write_file_neg, k, w, lamada)

    dataset_list = []
    for filename in file_list:
        temp = []
        with open(filename, 'r') as f:
            for line in f:
                temp.append(line.strip())
        dataset_list.append(temp)

    c_cost = 2 ** c
    g_gamma = 2 ** g
    svm_params += (' -c ' + str(c_cost) + ' -g ' + str(g_gamma))

    data_list = data_seperate(dataset_list, write_path, 5)
    merge_result = data_merge(data_list, write_path)
    average_rate = svm_process(merge_result[0], merge_result[1], metric, svm_params, bi_or_multi)
    lock.acquire()
    try:
        if average_rate > best_rate.value:
            best_rate.value = average_rate
            del best_params[:]
            best_params.extend([c, g])
            del t_v_list[:]
            t_v_list.extend(data_list)
        print 'Iteration', ' c = ', c, ' g = ', g, ' finished.'
        #print average_rate
    finally:
        lock.release()
    semph.release()


def data_seperate(dataset_list, path, k):
    """Seperate the datasets into k parts.
    :param pos: positive dataset, a list object.
    :param neg: negative dataset, a list object.
    :param path: the path of the generated files.
    :param k: int, the number of parts.
    """
    dataset_lens = []
    for dataset in dataset_list:
        dataset_lens.append(len(dataset))

    #len_pos = len(pos)
    #len_neg = len(neg)
    if min(dataset_lens) < k:
        k = min(dataset_lens)
    #if (len_pos < k) or (len_neg < k):
    #    k = min(len_pos, len_neg)

    size_of_part = [] # the size of every part of the k parts.
    for length in dataset_lens:
        size = (length + 1) / k
        size_of_part.append(size)

    sep_dataset_list = []
    for i in xrange(len(dataset_list)):
        rand_list = random_select(dataset_list[i], size_of_part[i], k)
        sep_dataset_list.append(rand_list)



    #len_part_pos = (len_pos + 1) / k
    #len_part_neg = (len_neg + 1) / k
    #ran_list_pos = random_select(pos, len_part_pos, k)
    #ran_list_neg = random_select(neg, len_part_neg, k)
    final_list = []


    while(len(sep_dataset_list[0]) != 0):
        temp_list = []
        for sep_dataset in sep_dataset_list:
            index = random.sample(range(len(sep_dataset)), 1)
            elem = sep_dataset.pop(index[0])
            temp_list.extend(elem)
        final_list.append(temp_list)



        # index_pos = random.sample(range(len(ran_list_pos)), 1)
        # index_neg = random.sample(range(len(ran_list_neg)), 1)
        # elem_pos = ran_list_pos.pop(index_pos[0])
        # elem_neg = ran_list_neg.pop(index_neg[0])
        # temp_list.extend(elem_pos)
        # temp_list.extend(elem_neg)
        # final_list.append(temp_list)
    #current_path = os.path.dirname(os.path.realpath(__file__))
    for count, part in enumerate(final_list):
        if count == 0:
            file_name = path + 'test_set.txt'
            with open(file_name, 'w') as f:
                for elem in part:
                    f.write(elem)
                    f.write('\n')
        else:
            file_name = path + 'seperate_set_' + str(count) + '.txt'
            with open(file_name, 'w') as f:
                for elem in part:
                    f.write(elem)
                    f.write('\n')

    return final_list


def random_select(data, size_of_part, k):
    """Divide dataset randomly.
    :param data: dataset to be seperated.
    :param size_of_part: size of every part.
    """
    rand_list = []
    count = 1
    while (count <= k - 1):
        part_list = []
        for i in range(size_of_part):
            index = random.sample(range(len(data)), 1)
            elem = data.pop(index[0])
            part_list.append(elem)
        rand_list.append(part_list)
        count += 1
    if len(data) != 0:
        rand_list.append(data)
    return rand_list



def data_merge(data_list, path):
    """merge seperated sets for trainning set.
    :param data_list: a list of the seperated dataset.
    :param path: the path of the generated files.
    """
    if len(data_list) != 0:
        test_list = data_list.pop(0)
    num_list = range(len(data_list))
    comb_list = list(combinations(num_list, len(data_list) - 1))
    training_list = []
    validation_list = []
    for elem in comb_list:
        train_set = []
        #vali_set = []
        for i in elem:
            train_set.extend(data_list[i])
        diff_num = list(set(num_list) - set(elem))
        validation_list.append(data_list[diff_num[0]])
        training_list.append(train_set)
    train_file_list = []
    #current_path = os.path.dirname(os.path.realpath(__file__))
    for count, part in enumerate(training_list):
        file_name = path + 'training_set_' + str(count + 1) + '.txt'
        with open(file_name, 'w') as f:
            for li in part:
                f.write(li)
                f.write('\n')
        train_file_list.append(file_name)
    validate_file_list = []
    for count, part in enumerate(validation_list):
        file_name = path + 'validation_set_' + str(count + 1) + '.txt'
        with open(file_name, 'w') as f:
            for li in part:
                f.write(li)
                f.write('\n')
        validate_file_list.append(file_name)

    return train_file_list, validate_file_list, training_list, validation_list, test_list



def svm_process(train_file_list, validate_file_list, metric, svm_params, bi_or_multi):
    """ The cross validation process for parameter selection.
    :param train_file_list: list of the training file names.
    :param validate_file_list: list of the validation file names.
    :param metric: the performance measurement used for parameter selection.
    :param svm_params: the parameters of svm_train in libsvm.
    """
    if len(train_file_list) != len(validate_file_list):
        raise ValueError("The number of the training files must equal to that of the validation files.")
    #t_label_list = []
    #p_label_list = []
    rate_list = []
    for i in range(len(train_file_list)):
        rates = one_svm_process(train_file_list[i], validate_file_list[i], metric, svm_params, bi_or_multi)
        if bi_or_multi == 0:
            rate_list.append(rates[0])
        elif bi_or_multi == 1:
            rate_list.append(rates)
    average = sum(rate_list) / len(rate_list)
    #for test
    #---------------------------------------
    #print rate_list, '\n'
    #---------------------------------------

    return average



def one_svm_process(train_file, validation_file, metric, svm_params, bi_or_multi):
    """This is the process of train & predict for one time.
    :param train_file: the training set file.
    :param validation_file: the validation set file.
    :param metric: the performance measurement used for parameter selection.
    :param svm_params: the parameters of svm_train in libsvm.
    """
    train_y, train_x = svm_read_problem(train_file)
    validation_label, prob_validation_data = svm_read_problem(validation_file)
    #svm_param = svm_parameter('-c 128 -g 0.5 -q')
    #svm_param = svm_parameter(svm_params)
    model = svm_train(train_y, train_x, svm_params)
    p_label, p_acc, p_val = svm_predict(validation_label, prob_validation_data, model, '-q')
    labels = model.get_labels()
    deci = [labels[0]*val[0] for val in p_val]
    evals = performance(validation_label, p_label, deci, bi_or_multi=bi_or_multi)
    #for test
    #---------------------------------------
    #print validation_label, '\n', p_label
    #---------------------------------------
    if bi_or_multi == 0:
        return evals[metric], evals
    elif bi_or_multi == 1:
        return evals




def performance(origin_labels, predict_labels, deci_value, output=None, title=None, roc=False, bi_or_multi=0):
    """evaluations used to evaluate the performance of the model.
    :param origin_labels: true values of the dataset.
    :param predict_labels: predicted values of the dataset.
    :param deci_value: decision values used for ROC and AUC.
    :param output: the output file name of the ROC curve.
    :param title: the title of the ROC curve.
    :param roc: indicate whether to draw the ROC curve or not.
    """
    if len(origin_labels) != len(predict_labels):
        raise ValueError("The number of the original labels must equal to that of the predicted labels.")
    if bi_or_multi == 0:
        TP = 0.0
        TN = 0.0
        FP = 0.0
        FN = 0.0
        for i in range(len(origin_labels)):
            if origin_labels[i] == 1.0 and predict_labels[i] == 1.0:
                TP += 1.0
            elif origin_labels[i] == 1.0 and predict_labels[i] == -1.0:
                FN += 1.0
            elif origin_labels[i] == -1.0 and predict_labels[i] == 1.0:
                FP += 1.0
            elif origin_labels[i] == -1.0 and predict_labels[i] == -1.0:
                TN += 1.0
        try:
            SN = TP / (TP + FN)
        except ZeroDivisionError:
            SN = 0.0
        try:
            SP = TN / (FP + TN)
        except ZeroDivisionError:
            SP = 0.0
        ACC = (TP + TN) / (TP + TN + FP + FN)
        try:
            MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        except ZeroDivisionError:
            MCC = 0.0
        current_path = os.path.dirname(os.path.realpath(__file__))
        roc_data_file = current_path + const.GEN_FILE_PATH + 'roc_data'
        try:
            AUC = plot_roc(deci_value, origin_labels, output, title, roc, roc_data_file)
        except ZeroDivisionError:
            AUC = 0.0
        del_file(roc_data_file)

        return ACC, MCC, AUC, SN, SP
    elif bi_or_multi == 1:
        correct_labels = 0.0
        for elem in zip(origin_labels, predict_labels):
            if elem[0] == elem[1]:
                correct_labels += 1.0
        ACC = correct_labels / len(origin_labels)
        return ACC


def cross_validation(label_list, vector_list, fold, svm_params, predict_params, bi_or_multi):
    """Do cross validation.
    :param label_list: list of labels.
    :param vector_list: list of vectors.
    :param fold: the fold of cross validation.
    """
    result = dataset_split_cv(label_list, vector_list, fold)
    if result == False:
        return False
    else:
        split_vector_list, split_label_list = result
    len_vector = len(split_vector_list)
    len_label = len(split_label_list)
    if len_vector != len_label:
        print 'Error: The length of the labels is not equal to that of the vectors.'
        return False
    deci = []
    acc_list = []
    mcc_list = []
    auc_list = []
    sn_list = []
    sp_list = []
    if bi_or_multi == 0:
        for i in range(len_vector):
            train_vector_list = []
            train_label_list = []
            #test_vector_list = []
            #test_label_list = []
            test_vector_list = split_vector_list[i]
            test_label_list = split_label_list[i]
            for j in range(len_vector):
                if j != i:
                    train_vector_list.extend(split_vector_list[j])
                    train_label_list.extend(split_label_list[j])
            m = svm_train(train_label_list, train_vector_list, svm_params)
            p_label, p_acc, p_val = svm_predict(test_label_list, test_vector_list, m, predict_params)
            labels = m.get_labels()
            subdeci = [labels[0]*val[0] for val in p_val]
            deci += subdeci
            evals = performance(test_label_list, p_label, subdeci, bi_or_multi=bi_or_multi)
            acc_list.append(evals[0])
            mcc_list.append(evals[1])
            auc_list.append(evals[2])
            sn_list.append(evals[3])
            sp_list.append(evals[4])
        acc_average = sum(acc_list) / len(acc_list)
        mcc_average = sum(mcc_list) / len(mcc_list)
        auc_average = sum(auc_list) / len(auc_list)
        sn_average = sum(sn_list) / len(sn_list)
        sp_average = sum(sp_list) / len(sp_list)

        label_all = []
        for i in split_label_list:
            label_all.extend(i)
        check_gnuplot_exe()
        roc_output = 'cross_validation.png'
        title = 'cross validation'
        current_path = os.path.dirname(os.path.realpath(__file__))
        roc_data_file = current_path + const.GEN_FILE_PATH + 'roc_data'
        plot_roc(deci, label_all, roc_output, title, True, roc_data_file)
        del_file(roc_data_file)
        dest_file = current_path + const.FINAL_RESULTS_PATH + roc_output
        copy_file(roc_output, dest_file)
        print ('The cross validation results are as follows:')
        print 'ACC = %.4f' % acc_average
        print 'MCC = %.4f' % mcc_average
        print 'AUC = %.4f' % auc_average
        print 'Sn  = %.4f' % sn_average
        print 'Sp  = %.4f\n' % sp_average
        print "The ROC curve has been saved. You can check it here: "
        if sys.platform.startswith('win'):
            print dest_file.replace('/', '\\'), '\n'
        else:
            print dest_file.replace('\\', '/'), '\n'
    elif bi_or_multi == 1:
        for i in range(len_vector):
            train_vector_list = []
            train_label_list = []
            #test_vector_list = []
            #test_label_list = []
            test_vector_list = split_vector_list[i]
            test_label_list = split_label_list[i]
            for j in range(len_vector):
                if j != i:
                    train_vector_list.extend(split_vector_list[j])
                    train_label_list.extend(split_label_list[j])
            m = svm_train(train_label_list, train_vector_list, svm_params)
            p_label, p_acc, p_val = svm_predict(test_label_list, test_vector_list, m, predict_params)
            labels = m.get_labels()
            subdeci = [labels[0]*val[0] for val in p_val]
            deci += subdeci
            evals = performance(test_label_list, p_label, subdeci, bi_or_multi=bi_or_multi)
            acc_list.append(evals)
        acc_average = sum(acc_list) / len(acc_list)
        print ('The cross validation results are as follows:')
        print 'ACC = %.4f' % acc_average




def dataset_split_cv(label_list, vector_list, fold):
    """Split dataset for cross validation.
    :param label_list: list of labels.
    :param vector_list: list of vectors.
    :param fold: the fold of cross validation.
    """
    length = len(label_list)
    if fold > length or fold <= 1:
        print 'Error: The fold of cross validation should be larger than 1 and less than or equal to the amount of vectors.'
        return False
    len_part = (length + 1) / fold
    split_vector_list = []
    split_label_list = []
    count = 1
    #index_list = range(length)
    while (count <= fold - 1):
        part_vector_list = []
        part_label_list = []
        for i in range(len_part):
            index = random.sample(range(len(vector_list)), 1)
            vector_elem = vector_list.pop(index[0])
            part_vector_list.append(vector_elem)
            label_elem = label_list.pop(index[0])
            part_label_list.append(label_elem)
        split_vector_list.append(part_vector_list)
        split_label_list.append(part_label_list)
        count += 1

    if len(vector_list) != 0:
        split_label_list.append(label_list)
        split_vector_list.append(vector_list)

    return split_vector_list, split_label_list


def check_args(args):
    """check the arguments of the command line.
    :param args: an object of the arguments.
    """
    if 'v' in args:
        if args.v is not None and args.v.isdigit() == True and int(args.v) < 2:
            print 'Error: If the value of -v is a digit, then it must be larger than 1.'
            return False
    return True

def check_c_g(arg, c_or_g):
    if arg is None:
        if c_or_g == 'c':
            x = 7
        elif c_or_g == 'g':
            x = -1
    elif len(arg) == 1:
        try:
            x = int(arg[0])
        except ValueError:
            print 'Error: the value of argument "-c" or "-g" should be a digit.'
            return False
    elif len(arg) == 2:
        try:
            a = int(arg[0])
            b = int(arg[1])
            x = xrange(int(min(a, b)), int(max(a, b))+1)
        except ValueError:
            print 'Error: the values of argument "-c" or "-g" should be digits.'
            return False
    elif len(arg) == 3:
        try:
            a = int(arg[0])
            b = int(arg[1])
            c = int(arg[2])
            x = xrange(int(min(a, b)), int(max(a, b))+1, int(c))
        except ValueError:
            print 'Error: the values of argument "-c" or "-g" should be digits.'
            return False
    elif len(arg) > 3:
        print 'Error: the number of values of "-c" or "-g" cannot be more than 3.'
        return False

    return x



def main(args):
    """The main process of train.
    :param args: an object of the arguments.
    """
    #================================================================================================
    #                   Check whether the command line arguments are valid or not.
    #================================================================================================
    start_time = time.time()
    # Path to find gnuplot.
    gnuplot_exe_list = [r'"C:\Program Files\gnuplot\pgnuplot.exe"',
                        r'".\gnuplot\bin\pgnuplot.exe"', "/usr/bin/gnuplot","/usr/local/bin/gnuplot"]
    # Get the current path.
    current_path = os.path.dirname(os.path.realpath(__file__))

    # Judge whether the path contains Chinese character or not.
    current_path_uni = unicode(current_path, "gbk")
    if check_contain_chinese(current_path_uni):
        print 'Error: the path can not contain Chinese characters.'
        return False

    file_list = args.files

    # Judge whether binary classification or multiclass classification.
    if len(file_list) == 2:
        bi_or_multi = 0
    elif len(file_list) > 2:
        bi_or_multi = 1
    else:
        print 'The number of input files must be more than 1.'
        return False

    preprocess_result = []
    for i in file_list:
        result = data_preprocess(i, const.TEMP_DIR)
        preprocess_result.append(result)

    if False in preprocess_result:
        print 'There exist some files that do not satisfy the LIBSVM format.'
        return False
    else:
        new_file_list = preprocess_result

    trans_labels(new_file_list)

    if args.v == 'i' and args.i_files is not None:
        for i in args.i_files:
            result = data_preprocess(i, const.TEMP_INDEPENDENT_DIR)
            preprocess_result.append(result)
        if False in preprocess_result:
            print 'There exist some independent test files that do not satisfy the LIBSVM format.'
            return False
        else:
            independent_file_list = preprocess_result


    predict_params = '-q' # optional parameters of svm_predict()

    svm_params = '-h 0 -m 1024 -q'

    c_result = check_c_g(args.c, 'c')
    if c_result is False:
        return False
    # elif type(c_result) == int:
    #     c = c_result
    # elif type(c_result) == xrange:
    #     c_range = c_result

    g_result = check_c_g(args.g, 'g')
    if g_result is False:
        return False
    # elif type(g_result) == int:
    #     g = g_result
    # elif type(g_result) == xrange:
    #     g_range = g_result
    if type(c_result) != type(g_result):
        print 'Both the arguments c and g should be specified values or both of them are ranges.'
        return False

    if args.b == '1':
        svm_params += (' -b ' + str(args.b))
        predict_params += (' -b ' + str(args.b))
        b = args.b
    elif args.b == '0':
        b = args.b

    if args.p == 'ACC':
        metric = 0
    elif args.p == 'MCC':
        metric = 1
    elif args.p == 'AUC':
        metric = 2
    if args.m is not None:
        model_file_name = args.m
    else:
        print 'Error: the name of the model can not be omitted.'
        print 'A value should be given to the parameter -m.'
        return False

    cpu_core = mul.cpu_count()
    if args.cpu is None:
        process_num = cpu_core
    elif 0 < args.cpu <= cpu_core:
        process_num = args.cpu
    elif args.cpu < 0 or args.cpu > cpu_core:
        process_num = cpu_core
        print 'Warning: The value of -cpu should be larger than 0'
        print 'and less than or equal to the number of cpu core in your computer.'
        print 'The value has been set as the default(number of all cpu cores in your computer).'
        time.sleep(2)

    if args.v == 'i' and args.i_files is None:
        print 'At least one independent dataset file should be included.'
        return False



    #================================================================================================
    #                                     Args check finished here.
    #================================================================================================

    #================================================================================================
    #                                    Parameter selection starts.
    #================================================================================================

    if type(c_result) == xrange:
        c_range = c_result
        g_range = g_result
        #print g_result

        print 'Parameter selection is in processing...\n'

        results = param_selection(new_file_list, metric, svm_params, process_num, c_range, g_range, bi_or_multi)

        print 'Parameter selection completed.\n'

        c = results[0][0]
        g = results[0][1]
        print 'The optimal parameters for the dataset are: c = ', c, ' g = ', g
        print '\n'

        if args.v is None:
            print 'The performance evaluations for the optimal parameter(s) are as follows:\n'
            if bi_or_multi == 0:
                print 'ACC = %.4f' % results[1][1][0]
                print 'MCC = %.4f' % results[1][1][1]
                print 'AUC = %.4f' % results[1][1][2]
                print 'Sn  = %.4f' % results[1][1][3]
                print 'Sp  = %.4f\n' % results[1][1][4]
            elif bi_or_multi == 1:
                print 'ACC = %.4f' % results[1]

    #================================================================================================
    #                                    Parameter selection finished.
    #================================================================================================

    elif type(c_result) == int:
        c = c_result
        g = g_result

    c_cost = 2 ** c
    g_gamma = 2 ** g
    svm_params += (' -c ' + str(c_cost) + ' -g ' + str(g_gamma))

    y_all = []
    x_all = []
    for file in new_file_list:
        y, x = svm_read_problem(file)
        y_all.extend(y)
        x_all.extend(x)

    dataset_size = len(x_all)
    pkl_y = current_path + const.TEMP_DIR.lstrip('.') + 'dataset_y.pkl'
    pkl_x = current_path + const.TEMP_DIR.lstrip('.') + 'dataset_x.pkl'
    cPickle.dump(y_all, open(pkl_y, 'wb'))
    cPickle.dump(x_all, open(pkl_x, 'wb'))


    #================================================================================================
    #                                    Model training & cross validation.
    #================================================================================================
    print 'Model training is in processing...'

    final_model_file = current_path + const.FINAL_RESULTS_PATH + model_file_name
    middle_model_file = current_path + const.FINAL_RESULTS_PATH + 'middle.model'
    # print y_all
    # print x_all
    # jackknife cross validation.
    if args.v =='j':
        cross_validation(y_all, x_all, dataset_size, svm_params, predict_params, bi_or_multi)
    # k-fold cross validation.
    elif args.v is not None and args.v.isdigit() == True and int(args.v) > 1:
        fold = int(args.v)
        cross_validation(y_all, x_all, fold, svm_params, predict_params, bi_or_multi)

    y_all = cPickle.load(open(pkl_y, 'rb'))
    x_all = cPickle.load(open(pkl_x, 'rb'))

    final_model = svm_train(y_all, x_all, svm_params)
    svm_save_model(middle_model_file, final_model)

    #================================================================================================
    #                             Add the parameters to the SVM model file.
    #================================================================================================

    middle_list = []
    with open(middle_model_file) as f:
        for i in f:
            middle_list.append(i)

    param_line = 'c,' + str(c) + ',g,' + str(g) + ',b,' + str(b) + ',bi_or_multi,' + str(bi_or_multi)
    with open(final_model_file, 'w') as f:
        f.write(param_line)
        f.write('\n')
        for i in middle_list:
            f.write(i)

    print 'Model training completed.'
    print 'The model has been saved. You can check it here:'
    if sys.platform.startswith('win'):
        print final_model_file.replace('/', '\\'), '\n'
    else:
        print final_model_file.replace('\\', '/'), '\n'

    if os.path.isfile('cross_validation.png'):
        try:
            os.remove('cross_validation.png')
        except OSError:
            time.sleep(0.1)
            try:
                os.remove('cross_validation.png')
            except OSError:
                pass

    #================================================================================================
    #                                   Independent dataset test.
    #================================================================================================

    if 'independent_file_list' in locals().keys():
        print 'The independent test dataset is found.\n'
        test_y = []
        test_x = []
        for file in independent_file_list:
            y, x = svm_read_problem(file)
            test_y.extend(y)
            test_x.extend(x)
        model = svm_load_model(middle_model_file)
        p_label, p_acc, p_val = svm_predict(test_y, test_x, model, predict_params)
        labels = model.get_labels()
        deci = [labels[0]*val[0] for val in p_val]
        check_gnuplot_exe()
        roc_output = 'independent_roc.png'
        title = 'the test dataset'
        evals = performance(test_y, p_label, deci, roc_output, title, True, bi_or_multi)
        if bi_or_multi == 0:
            dest_file = current_path + const.FINAL_RESULTS_PATH + roc_output
            copy_file(roc_output, dest_file)
            print 'The performance evaluations of the final model are as follows:\n'
            print 'ACC = %.4f' % evals[0]
            print 'MCC = %.4f' % evals[1]
            print 'AUC = %.4f' % evals[2]
            print 'Sn  = %.4f' % evals[3]
            print 'Sp  = %.4f\n' % evals[4]
            print "The ROC curve has been saved. You can check it here: "
            if sys.platform.startswith('win'):
                print dest_file.replace('/', '\\'), '\n'
            else:
                print dest_file.replace('\\', '/'), '\n'

            if os.path.isfile('independent_roc.png'):
                try:
                    os.remove('independent_roc.png')
                except OSError:
                    time.sleep(0.1)
                    try:
                        os.remove('independent_roc.png')
                    except OSError:
                        pass
        elif bi_or_multi == 1:
            print 'The performance evaluations of the final model are as follows:\n'
            print 'ACC = %.4f' % evals
    print '\n'
    print 'Done.'
    print 'Used time: %.2fs' % (time.time() - start_time)


if __name__ == '__main__':
    #data_preprocess('res1.txt')
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="The train module for training a SVM classifier.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('files', nargs='*',
                       help="The input files, generated by former steps or webserver.")
    parse.add_argument('-p', default='ACC', choices=['ACC', 'MCC', 'AUC'],
                       help="The performance metric used for parameter selection.\n"
                       "Default value is ACC.")
    parse.add_argument('-m',
                       help="The name of the trained model.")
    parse.add_argument('-cpu', type=int,
                       help="The maximum number of CPU core used for multiprocessing in\n"
                       "parameter selection. Default is the number of all CPU cores.")
    parse.add_argument('-c', nargs='*',
                       help="The parameter cost of RBF kernel is usually represented as the\n"
                       " c power of 2 and the value of c should be input here. Either one value \n"
                       "or a range is acceptable. If users want to input a range, the format is like\n"
                       " '-c 1 5 1': the first value is the lower bound, the second value is the upper\n"
                       " bound and the third value is the step.\n"
                       "Default value is 7.")
    parse.add_argument('-g', nargs='*',
                       help="The parameter gamma of RBF kernel is usually represented as the\n"
                       " g power of 2 and the value of g should be input here. Either one value \n"
                       "or a range is acceptable. If users want to input a range, the format is like\n"
                       " '-g -1 5': the first value is the lower bound, the second value is the upper\n"
                       " bound and the third value is the step. The default step is 1.\n"
                       "Default value is -1.")
    parse.add_argument('-b', default='0', choices=['0', '1'],
                       help="whether to train a SVC or SVR model for probability\n"
                       "estimates, 0 or 1. Default value is 0.")
    parse.add_argument('-v',
                       help="The cross validation mode.\n"
                       "n: (an integer larger than 0) n-fold cross validation.\n"
                       "j: (character 'j') jackknife cross validation.\n"
                       "i: (character 'i') independent test set method.")
    parse.add_argument('-i_files', nargs='+',
                       help="The independent test dataset.\n"
                       "If the parameter '-v' is specified as 'i', one or more\n"
                       " independent test dataset files should be\n"
                       "included.\n"
                       "e.g. '-i_files test1.txt test2.txt'.")
    args = parse.parse_args()
    if check_args(args):
        print 'Processing...'
        main(args)
