import pickle

## Reads and prints the global_embeddings dataset

def main():
    with open('Datasets/global_embeddings_map.pkl', 'rb') as f:
        data = pickle.load(f)
        print(data)

if __name__ == "__main__":
    main()