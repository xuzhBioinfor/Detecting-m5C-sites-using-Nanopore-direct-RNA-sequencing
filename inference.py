import os, sys
import torch
import toml
import warnings
import numpy as np
from model.model import MILModel
# from constants import DEFAULT_MODEL_CONFIG, DEFAULT_MODEL_WEIGHTS, DEFAULT_NORM_PATH, DEFAULT_MIN_READS, DEFAULT_READ_THRESHOLD
from data_utils import NanopolishDS, NanopolishReplicateDS, inference_collate
from inference_utils import run_inference
from torch.utils.data import DataLoader
import argparse
from datetime import datetime


def run(input_dir, out_dir, model_config, device, model_state_dict, norm_path, n_processes, batch_size, save_per_batch,
        read_proba_threshold, num_iterations,DEFAULT_MIN_READS):
    model = MILModel(toml.load(model_config)).to(device)
    model.load_state_dict(torch.load(model_state_dict, map_location=torch.device(device)))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(os.path.join(out_dir, "data.site_proba.csv"), 'w', encoding='utf-8') as f:
        f.write('transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio\n')
    with open(os.path.join(out_dir, "data.indiv_proba.csv"), 'w', encoding='utf-8') as g:
        g.write('transcript_id,transcript_position,read_index,probability_modified\n')

    if len(input_dir) == 1:
        ds = NanopolishDS(input_dir[0], DEFAULT_MIN_READS, norm_path, mode='Inference')
    else:
        ds = NanopolishReplicateDS(input_dir, DEFAULT_MIN_READS, norm_path, mode='Inference')

    dl = DataLoader(ds, num_workers=n_processes, collate_fn=inference_collate, batch_size=batch_size, shuffle=False)
    # run_inference(model, dl, args)
    run_inference(model, dl, out_dir, device, save_per_batch, read_proba_threshold, num_iterations, n_processes)


def proprecess_input(inputs, outputs, threads, seed, device, batch_size, save_per_batch, num_iterations,
                     DEFAULT_MIN_READS):
    # input dir
    if os.path.lexists(os.path.expanduser(inputs)) is False:
        sys.exit("[ERROR]Can't get your input dir, please check the PATH!\n")
    else:
        inputs = [inputs]
    sys.stderr.write("Input dir is ready!\n")

    # output
    if os.path.lexists(os.path.expanduser(outputs)) is True and \
            os.path.isdir(os.path.expanduser(outputs)) is True:
        sys.stderr.write("Output path is ready!\n")
        print("OUTPUT PATH : %s" % (outputs))
    else:
        outputs = os.getcwd()
        print("[WARNING]Can't get your save path, the result file will be output in the WORK PATH!\n")
        print("WORK PATH : %s" % (outputs))

    # threads
    if isinstance(threads, int):
        print("Threads : %d" % threads)
    else:
        sys.stderr.write("Threads is not an INTEGER!\n")

    # seed
    if isinstance(seed, int):
        print("Seed : %d" % seed)
    else:
        sys.stderr.write("Seed is not an INTEGER!\n")

    # device (cpu/gpu)
    device = device.lower()
    if device == 'cpu' or device == 'gpu':
        print("Device : %s" % device)
    else:
        device = 'cpu'
        print("[WARNING]Device not in ‘cpu’ and 'gpu'! Device will be set to 'cpu'！")
        print("Device : %s" % device)

    # batch_size
    if isinstance(batch_size, int):
        print("Batch size : %d" % batch_size)
    else:
        sys.stderr.write("Batch size is not an INTEGER!\n")

    # save_per_batch
    if isinstance(save_per_batch, int):
        print("Save per batch : %d" % save_per_batch)
    else:
        sys.stderr.write("Save per batch is not an INTEGER!\n")

    # num_iterations
    if isinstance(num_iterations, int):
        print("Num iterations : %d" % num_iterations)
    else:
        sys.stderr.write("Num iterations is not an INTEGER!\n")

    # DEFAULT_MIN_READS
    if isinstance(DEFAULT_MIN_READS, int):
        print("DEFAULT MIN READS : %d" % DEFAULT_MIN_READS)
    else:
        sys.stderr.write("DEFAULT MIN READS is not an INTEGER!\n")

    # default settings
    model_config = os.getcwd() + '/model/configs/model_configs/m5C_mouse_model.toml'  # .toml
    model_state_dict = os.getcwd() +'/model/model_states/m5C_mouse_model.pt'  # .pt
    norm_path = os.getcwd() +'/model/norm_factors/m5C_mouse_norm_dict_nanopolish.joblib'  # .joblib
    read_proba_threshold = 0.1634889841079712

    print(inputs)

    run(inputs, outputs, model_config, device, model_state_dict, norm_path, threads, batch_size, save_per_batch,
        read_proba_threshold, num_iterations,DEFAULT_MIN_READS)


def main():
    parser = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input_dir", required=1, type=str, help="input file path")
    parser.add_argument("-o", "--out_dir", required=1, type=str, help="output path")
    parser.add_argument("-t", "--threads", default=1, type=int, help="number of cores to use [default:1]")
    parser.add_argument("-s", "--seed", default=0, type=int, help="[default:0]")
    parser.add_argument("-d", "--device", default='cpu', type=str, help="cpu or gpu[default:cpu]")
    parser.add_argument("-b", "--batch_size", default=512, type=int, help="[default:512]")
    parser.add_argument("-p", "--save_per_batch", default=2, type=int, help="[default:2]")
    parser.add_argument("-n", "--num_iterations", default=1000, type=int, help="[default:1000]")
    parser.add_argument("-m", "--DEFAULT_MIN_READS", default=20, type=int, help="[default:20]")

    o = parser.parse_args()
    proprecess_input(o.input_dir, o.out_dir, o.threads, o.seed, o.device, o.batch_size,
                     o.save_per_batch, o.num_iterations, o.DEFAULT_MIN_READS)


if __name__ == '__main__':
    desc = "Using the output result from 'dataprep.py' to infer the site-level probability and methylation ratio "
    epilog = "This step is modified from m6Anet. \
    You can see this project in github: https://github.com/xuzhBioinfor/Detecting-m5C-sites-using-Nanopore-direct-RNA-sequecing \
    and m6Anet in github: https://github.com/GoekeLab/m6anet"

    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!\n")

    dt = datetime.now() - t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
