import pickle

def pickle_object(obj, filename, protocol=pickle.HIGHEST_PROTOCOL):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, protocol = protocol)
    return

def unpickle_object(filename):
    with open(filename, 'rb') as f:
        obj = pickle.load(f)
    return obj