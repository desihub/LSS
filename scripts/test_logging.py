import logging

# create logger
logger = logging.getLogger('test')
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

# 'application' code
logger.debug('debug message')
logger.info('info message')
#logger.warning('warn message')
#logger.error('error message')
#logger.critical('critical message')

#import logging

#logger = logging.getLogger('test')
#logger.setLevel(level=logging.INFO)
#logger.info('info message')