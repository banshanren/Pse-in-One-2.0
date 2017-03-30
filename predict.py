#!/usr/bin/env python
# -*- coding: utf-8 -*-


import time, os, sys
from util import check_contain_chinese, copy_file
import const
from train import performance, data_preprocess
from libsvm.python.svmutil import *
from libsvm.python.plotroc import *



def main(args):
    """The main process of predict.
    :param args: an object of the arguments.
    """
    start_time = time.time()

    current_path = os.path.dirname(os.path.realpath(__file__))

    #Judge whether the path contains Chinese character or not
    current_path_uni = unicode(current_path, "gbk")
    if check_contain_chinese(current_path_uni):
        print 'Error: the path can not contain Chinese character.'
        return False

    #================================================================================================
    #                                      Inputfile preprocess.
    #================================================================================================

    inputfile = args.inputfile

    result = data_preprocess(inputfile, const.TEMP_DIR)
    if result == False:
        print 'The input file does not satisfy the LIBSVM format.'
        return False
    else:
        new_inputfile = result


    if args.m == None:
        print 'Error: the name of the model can not be omitted.'
        print 'A value should be given to the parameter -m.'
        return False
    else:
        model_name = args.m
    model_file = current_path + const.FINAL_RESULTS_PATH + model_name

    if args.o == None:
        output_name = 'output_labels.txt'
    else:
        output_name = args.o
    output = current_path + const.FINAL_RESULTS_PATH + output_name


    #================================================================================================
    #                   Processing the model file generated in the train step.
    #================================================================================================

    param_dict = dict()
    model_list = []
    with open(model_file) as f:
        train_params = f.readline().strip()
        for line in f:
            model_list.append(line)

    svm_model_file = current_path + const.FINAL_RESULTS_PATH + 'svm_model.model'
    with open(svm_model_file, 'w') as f:
        for i in model_list:
            f.write(i)

    param_list = train_params.split(',')
    for index in range(0, len(param_list), 2):
        param_dict[param_list[index]] = param_list[index+1]

    if 'c' in param_dict.keys() and 'g' in param_dict.keys() and 'b' in param_dict.keys() and 'bi_or_multi' in param_dict.keys():
        c = int(param_dict['c'])
        g = int(param_dict['g'])
        b = param_dict['b']
        bi_or_multi = int(param_dict['bi_or_multi'])

    print 'The parameters of RBF kernel:'
    print 'c = ', c, ' g = ', g


    #================================================================================================
    #                                       Predicting process.
    #================================================================================================

    label_list = []
    if args.labels !=None:
        with open(args.labels) as f:
            for i in f:
                if i.strip() == '+1':
                    label_list.append(1.0)
                elif i.strip() == '-1':
                    label_list.append(-1.0)
                else:
                    label_list.append(int(i.strip()))

    predict_params = '-q'
    if b == '1':
        predict_params += (' -b ' + b)

    y ,x = svm_read_problem(new_inputfile)
    #print y
    model = svm_load_model(svm_model_file)
    model_labels = model.get_labels()
    p_label, p_acc, p_val = svm_predict(y, x, model, predict_params)
    #print p_label
    if bi_or_multi == 0:
        with open(output, 'w') as f:
            for i in p_label:
                if i == 1.0:
                    f.write('+1')
                if i == -1.0:
                    f.write('-1')
                f.write('\n')
        if len(label_list) != 0:
            check_gnuplot_exe()
            deci = [model_labels[0]*val[0] for val in p_val]
            roc_output = 'predicted_roc.png'
            title = 'the predicted dataset'
            #print '1'
            evals = performance(label_list, p_label, deci, roc_output, title, True, bi_or_multi)
            dest_file = current_path + const.FINAL_RESULTS_PATH + roc_output
            copy_file(roc_output, dest_file)
            print 'The performance evaluations are as follows:\n'
            print 'ACC = %.4f' % evals[0]
            print 'MCC = %.4f' % evals[1]
            print 'AUC = %.4f' % evals[2]
            print 'Sn  = %.4f' % evals[3]
            print 'Sp  = %.4f\n' % evals[4]
            if evals[2] !=0:
                print "The ROC curve has been saved. You can check it here: "
                if sys.platform.startswith('win'):
                    print dest_file.replace('/', '\\'), '\n'
                else:
                    print dest_file.replace('\\', '/'), '\n'

        if os.path.isfile('predicted_roc.png'):
            try:
                os.remove('predicted_roc.png')
            except OSError:
                time.sleep(0.1)
                try:
                    os.remove('predicted_roc.png')
                except OSError:
                    pass

    elif bi_or_multi == 1:
        with open(output, 'w') as f:
            for i in p_label:
                f.write(str(i))
                f.write('\n')
        if len(label_list) != 0:
            deci = [model_labels[0]*val[0] for val in p_val]
            roc_output = 'predicted_roc.png'
            title = 'the predicted dataset'
            #print '1'
            evals = performance(label_list, p_label, deci, roc_output, title, True, bi_or_multi)
            print 'The performance evaluation is as follow:\n'
            print 'ACC = %.4f' % evals

    print "The predicted labels have been saved. You can check it here: "
    if sys.platform.startswith('win'):
        print output.replace('/', '\\'), '\n'
    else:
        print output.replace('\\', '/'), '\n'


    print("Done.")
    print("Used time: %.2fs" % (time.time() - start_time))





if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="The predict module for predicting a input file with a trained SVM classifier.",
                                    formatter_class=RawTextHelpFormatter)

    parse.add_argument('inputfile',
                       help="The input feature vectors file, in LIBSVM format.")
    parse.add_argument('-m',
                       help="The name of trained model generated by using 'train.py'.")
    parse.add_argument('-labels',
                       help="The real label file. Optional.")
    parse.add_argument('-o',
                       help="The output file name listing the predicted labels.\n"
                       "The default name is 'output_labels.txt'.")

    args = parse.parse_args()
    print 'Processing...'
    main(args)