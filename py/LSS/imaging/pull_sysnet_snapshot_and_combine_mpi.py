import sys
import os
from sysnet.sources.utils import split_NtoM
from sysnet.sources.train import forward
from sysnet.sources.models import init_model
from sysnet.sources.io import load_checkpoint, ImagingData, MyDataSet, DataLoader

from time import time
import fitsio as ft
import numpy as np
from astropy.table import Table
from glob import glob

from mpi4py import MPI

def chck2pid(chck):
    ch_ = chck.split('/')
    return '_'.join([ch_[-2], ch_[-1].split('.')[0]])

def do_forward(templates, checkpoints, rank, oudir, model_name, nnstruct, nfolds=5):

    axes = [i for i in range(templates['features'].shape[1])] # this will use all axes
    num_features = len(axes) 
    
    model = init_model(model_name)
    model = model(*nnstruct, input_dim=num_features)
    
    for i, chck in enumerate(checkpoints):        

        t0 = time()
        checkpoint = load_checkpoint(chck, model)
        img_data = ImagingData(templates, checkpoint['stats'], axes=axes)
        dataloader = DataLoader(MyDataSet(img_data.x, img_data.y, img_data.p, img_data.w),
                                 batch_size=2000000,
                                 shuffle=False) # num_workers=4
                                
        if rank == 0:print('finish data', time()-t0, i)
        result = forward(model, dataloader, {'device':'cpu'})        
        nnw = result[1].numpy().flatten()
        hpix = result[0].numpy()

        pid = chck2pid(chck)
        ouname = f'{oudir}/window_{pid}.fits'

        if rank == 0:print('finish forward pass ', time()-t0, i)
        
        dt = Table([hpix, nnw], names=['hpix', 'weight'])
        dt.write(ouname, format='fits')
        if rank == 0:print(f'save in {ouname}')
         
if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--type", help="tracer type to be selected")
    parser.add_argument('-i',  '--input_directory', type=str, help='path to the directory that contains the input templates')
    parser.add_argument('-sd', '--snapshot_directory',  type=str, help='path to the directory that contains the snapshots')
    parser.add_argument('-o',  '--output_directory', type=str, help='path to the parent directory to save windows')
    parser.add_argument('--north_nn_structure', type=int, nargs='*', default=(3, 10), help='North NN structure ( # hidden layer, # neurons)')
    parser.add_argument('--south_nn_structure', type=int, nargs='*', default=(4, 20), help='South NN structure ( # hidden layer, # neurons)')
    parser.add_argument('--model', type=str, default='dnnp', choices=['dnn','dnnp'], help='model, dnn or dnnp')
    args = parser.parse_args()
    
    zrl = [[0.8,1.1],[1.1,1.6]]
    regions = ['N','S']

    # if rank == 0:
    #     print(args)
    #     templates_dict = dict()
    #     chcks_dict = dict()
    #     oudir_dict = dict()
    #     for zl in zrl:
    #         zw = str(zl[0])+'_'+str(zl[1])
    #         for reg in regions:
    #             temp_id = f'{zw}_{reg}'
    #             templates_dict[temp_id] = ft.read(args.input_directory+f'/prep_{args.type}{zw}_{reg}.fits')
    #             chcks_dict[temp_id] = glob(args.snapshot_directory+f'/{args.type}{zw}_{reg}/model_*/snapshot_*.pth.tar')
    #             oudir_dict[temp_id] = args.output_directory+f'/windows_{args.type}{zw}_{reg}/'
    #             if not os.path.exists(oudir_dict[temp_id]):
    #                 os.makedirs(oudir_dict[temp_id])
    #                 print('rank 0 creates ', oudir_dict[temp_id])
    # else:
    #     templates_dict = None
    #     chcks_dict = None
    #     oudir_dict = None
    
    # templates_dict = comm.bcast(templates_dict, root=0)
    # chcks_dict = comm.bcast(chcks_dict, root=0)
    # oudir_dict = comm.bcast(oudir_dict, root=0)
    
    # for zl in zrl:
    #     zw = str(zl[0])+'_'+str(zl[1])
    #     for reg in regions:
    #         temp_id = f'{zw}_{reg}'
    #         chcks = chcks_dict[temp_id]
    #         templates = templates_dict[temp_id]
    #         oudir = oudir_dict[temp_id]
    #         nn_structure = args.north_nn_structure if reg == 'N' else args.south_nn_structure
            
    #         start, end = split_NtoM(len(chcks), size, rank)
    #         my_chcks = chcks[start:end+1]

    #         if rank == 0:print(rank, len(my_chcks), templates.size, my_chcks[:2])
    #         do_forward(templates, my_chcks, rank, oudir, args.model, nn_structure)
            
    # comm.Barrier()
    
    if rank == 0:
        # Generate and save nn-weights.fits per data split. 
        for zl in zrl:
            zw = str(zl[0])+'_'+str(zl[1])
            for reg in regions:
                temp_id = f'{args.type}{zw}_{reg}'
                windows_ = glob(args.output_directory+f'/windows_{temp_id}/window*')
        
                hpix = [] # hpix from first iteration
                weights_list = [] # list of weights 
                for i,iwindow in enumerate(windows_):
                    iw = ft.read(iwindow)
                    if i == 0: 
                        hpix = iw['hpix']
                    if not np.allclose(hpix, iw['hpix']):
                        raise ValueError('HEALPix pixels do not match!')
                    else:
                        weights_list.append(iw['weight'])
        
                shape = (len(hpix),len(windows_))
                print(shape)
                weights = np.zeros(shape[0], dtype=[('hpix', 'i8'), ('weight', 'f8', (shape[1], ))])
                weights['hpix']   = hpix
                weights['weight'] = np.stack(weights_list,axis=1)
                weights_path = args.output_directory + f'/windows_{temp_id}/nn-weights.fits'
                ft.write(weights_path, weights, clobber=True)
                print(f"wrote {weights_path}.")


    
    

    
    