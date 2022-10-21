import os
import argparse
import papermill     as pm

from   bin.pipeline  import pipeline
from   tidyup        import tidyup
from   findfile      import fetch_fields
from   pathlib       import Path


# https://docs.pytest.org/en/6.2.x/
def test_allnbs(survey='gama', no_qa=False):
    if 'GITHUB_ACTIONS' in os.environ:
        os.environ['USER']         = 'Hal' 
        os.environ['CODE_ROOT']    = os.environ['GITHUB_WORKSPACE']
        os.environ['GOLD_DIR']     = 'GAMA4/'
        os.environ['GOLD_LOGS']    = 'GAMA4/logs/'
        os.environ['RANDOMS_DIR']  = 'GAMA4/randoms/'
        os.environ['RANDOMS_LOGS'] = 'GAMA4/randoms/logs/'

        os.environ['PATH']         = os.environ['GITHUB_WORKSPACE'] + ':' + os.environ['GITHUB_WORKSPACE'] + '/bin:' + os.environ['PATH']
        os.environ['PYTHONPATH']   = os.environ['GITHUB_WORKSPACE'] + ':' + os.environ['GITHUB_WORKSPACE'] + '/bin:' + os.environ['PATH']

        Path(os.environ['GOLD_LOGS']).mkdir(parents=True, exist_ok=True)
        Path(os.environ['RANDOMS_LOGS']).mkdir(parents=True, exist_ok=True)

    if no_qa | ('GITHUB_ACTIONS' in os.environ):
        pipeline(use_sbatch=False, reset=False, nooverwrite=False, dryrun=True, survey='gama', freshclone=False, args=None)
        
        return 0

    print('Running all tests.')
        
    if (survey != 'gama') and (survey != 'desi'):
        raise  NotImplementedError(f'No implementation for survey: {survey}')
    
    tidyup()
    
    run_randomqa(survey)

    run_goldqa(survey)
    
    run_delta8qa(survey)
    
    print('Done.')

    return 0
    
def run_randomqa(survey):
    fields = fetch_fields(survey)

    print(fields)

    for field in fields:
        print('Running random QA for field {}'.format(field))

        pm.execute_notebook('docs/nb/randoms_n8_qa.ipynb',\
                            'test/pm_randoms_n8_{}_qa.ipynb'.format(field),\
                            parameters=dict(field=field,survey=survey),\
                            kernel='lumfn',\
                            )

def run_goldqa(survey):
    print('Running gold QA')

    tests = ['zmax_catQA', 'kE_catQA', 'ddp_QA', 'lumfn', 'delta8_qa']

    for test in tests:
        try:
            pm.execute_notebook('docs/nb/{}.ipynb'.format(test),\
                                'test/pm_{}.ipynb'.format(test),\
                                parameters=dict(survey=survey),\
                                kernel='lumfn',\
            )

        except:
            print('Failed on {} test with error: '.format(test))
            #print(E)

    '''
    pm.execute_notebook('docs/nb/kE_catQA.ipynb',\
                        'test/pm_kE_catQA.ipynb',\
                        kernel='lumfn',\
                        )

    pm.execute_notebook('docs/nb/ddp_QA.ipynb',\
                        'test/pm_ddp_QA.ipynb',\
                        kernel='lumfn',\
                        )

    pm.execute_notebook('docs/nb/lumfn.ipynb',\
                        'test/pm_lumfn.ipynb',\
                        kernel='lumfn',\
                        )

    pm.execute_notebook('docs/nb/delta8_qa.ipynb',\
                        'test/pm_delta8_qa.ipynb',\
                        kernel='lumfn',\
                        )
    '''
    # jack knife qa. 
    # d8 LF. 
    # desi qa. 
    
def run_delta8qa(survey):
    fields = fetch_fields(survey)
    
    for field in fields:
        print('Running delta8 QA for field {}'.format(field))
        
        pm.execute_notebook('docs/nb/d8LF_qa.ipynb',\
                    'test/pm_delta8_qa_{}.ipynb'.format(field),\
                    parameters=dict(field=field, survey=survey),\
                    kernel='lumfn',\
                    )

if __name__ == '__main__':
    # python3 -m pytest test_allnb.py 
    test_allnbs(no_qa=False)
