from    subprocess import check_output


def run_command(cmd, noid=False):
    print('Command: {}'.format(cmd))

    cmd = cmd.split()

    env = {}

    # env.update(os.environ)                                                                                                                                                                               
    # print('Calling ...')                                                                                                                                                                                  
    out = check_output(cmd)
    out = out.decode('utf-8')
    out = out.replace('\n', '')

    if noid:
        out=0

    out = int(out)

    return out
