import argparse
import os
import requests
def parseRescue(rescuefile):
    #Rescuefile parsing rules
    #Total number of node Total number of Nodes: XXXX
    #Number of failed node is written in format of Nodes that failed: XXXX
    #and below that line, name of failed node in written as # failed node 1, failed node 2, failed node 3, ....
    nTotalNodes = 0
    nFailedNodes = 0
    failedNodes = []
    
    with open(rescuefile, 'r') as f:
        lines = f.readlines()
    failedNodeInfo = ""
    for i, line in enumerate(lines):
        if 'Total number of Nodes' in line:
            nTotalNodes = int(line.split(':')[-1].strip())
        if 'Nodes that failed' in line:
            nFailedNodes = int(line.split(':')[-1].strip())
            failedNodeInfo = lines[i+1].replace('#', '')
            break
    failedNodes = failedNodeInfo.split(',')
    
    return nTotalNodes, nFailedNodes, failedNodes
    
if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--TOKEN', type=str, required=True)
    argparse.add_argument('--chatID', type=str, required=True)
    argparse.add_argument('--master_dir', type=str, required=True)
    args = argparse.parse_args()
    
    
    
    master_dir = args.master_dir
    sendingMessage = f'''Your Job is finished in the following directory: {master_dir} is finished\n'''
    sendingMessage += f'''Check All Jobs are finished without error....\n'''

    dag_dir = os.path.join(master_dir, 'dags')
    filelist = os.listdir(dag_dir)
    #check rescue file
    if not any('dagfile.dag.rescue' in file for file in filelist):
        sendingMessage += f'''No rescue file is generated by DAGMan\n'''
        
    else:
        sendingMessage += f'''Rescue file is generated by DAGMan\n'''
        rescue_file = [file for file in filelist if 'dagfile.dag.rescue' in file][-1]
        nTotalNodes, nFailedNodes, failedNodes = parseRescue(os.path.join(dag_dir, rescue_file))
        sendingMessage += f'''Total number of nodes: {nTotalNodes}\n'''
        sendingMessage += f'''Number of failed nodes: {nFailedNodes}\n'''
        for i,failedNode in enumerate(failedNodes):
            if 'ENDLIST' in failedNodes:
                    break
            sendingMessage += f'''Failed node: {failedNode}\n'''
            if i == 10:
                
                    
                sendingMessage += f'''And more....\n'''
                break
            
    requests.get(f'https://api.telegram.org/bot{args.TOKEN}/sendMessage?chat_id={args.chatID}&text={sendingMessage}')
        