import numpy as np
from keras.layers import Input, Dense, Activation, BatchNormalization, Dropout
from keras.models import Model
from keras.utils import to_categorical
from sklearn.datasets import make_blobs
from sklearn.metrics import accuracy_score
from keras import  optimizers
import matplotlib.pyplot as plt


def build_models(inp_dim, emb_dim, n_cls_source, alpha=2, alpha_lr=10):
    inputs = Input(shape=(inp_dim,)) 
    x4 = Dense(256, activation='linear')(inputs)
    x4 = BatchNormalization()(x4)
    x4 = Activation("elu")(x4)  
    x4 = Dense(emb_dim, activation='linear')(x4)
    x4 = BatchNormalization()(x4)
    x4 = Activation("elu")(x4)  

    source_classifier = Dense(n_cls_source, activation='softmax', name="mo1")(x4)  
    source_classifier = Activation('softmax', name='mo')(source_classifier)
    domain_classifier = Dense(32, activation='linear', name="do4")(x4)
    domain_classifier = BatchNormalization(name="do5")(domain_classifier)
    domain_classifier = Activation("elu", name="do6")(domain_classifier)
    domain_classifier = Dropout(0.5)(domain_classifier)

    domain_classifier = Dense(2, activation='softmax', name="do")(domain_classifier)

    comb_model = Model(inputs=inputs, outputs=[source_classifier, domain_classifier])
    comb_model.compile(optimizer="Adam",
              loss={'mo': 'categorical_crossentropy', 'do': 'categorical_crossentropy'},
              loss_weights={'mo': 1, 'do': alpha}, metrics=['accuracy'], )

    source_classification_model = Model(inputs=inputs, outputs=[source_classifier])
    source_classification_model.compile(optimizer=optimizers.adam(lr=0.001),
              loss={'mo': 'categorical_crossentropy'}, metrics=['accuracy'], )


    domain_classification_model = Model(inputs=inputs, outputs=[domain_classifier])
    domain_classification_model.compile(optimizer=optimizers.adam(lr=alpha_lr*0.001),
                  loss={'do': 'categorical_crossentropy'}, metrics=['accuracy'])
    
    
    embeddings_model = Model(inputs=inputs, outputs=[x4])
    embeddings_model.compile(optimizer="Adam",loss = 'categorical_crossentropy', metrics=['accuracy'])
                        
                        
    return comb_model, source_classification_model, domain_classification_model, embeddings_model

def batch_generator(data, batch_size):
    """Generate batches of data.

    Given a list of numpy data, it iterates over the list and returns batches of the same size
    This
    """
    all_examples_indices = len(data[0])
    while True:
        mini_batch_indices = np.random.choice(all_examples_indices, size=batch_size, replace=False)
        tbr = [k[mini_batch_indices] for k in data]
        yield tbr
        
        

def train( Xs, ys, Xt, yt=None, 
          emb_dim=2,
          batch_size = 64, 
          enable_dann = True, 
          n_iterations = 1000,
          alpha=2,
          alpha_lr=10):
    
    
    inp_dim = Xs.shape[1]
    ncls_source = len(set(ys))
    
    model, source_classification_model, domain_classification_model, embeddings_model = \
          build_models(inp_dim, emb_dim, ncls_source, alpha=alpha, alpha_lr = alpha_lr)

    y_adversarial_1 = to_categorical(np.array(([1] * batch_size + [0] * batch_size)))
    
    sample_weights_class = np.array(([1] * batch_size + [0] * batch_size))
    sample_weights_adversarial = np.ones((batch_size * 2,))

    S_batches = batch_generator([Xs, to_categorical(ys)], batch_size)
    T_batches = batch_generator([Xt, np.zeros(shape = (len(Xt),2))], batch_size)
    
    for i in range(n_iterations):
        # # print(y_class_dummy.shape, ys.shape)
        y_adversarial_2 = to_categorical(np.array(([0] * batch_size + [1] * batch_size)))

        X0, y0 = next(S_batches)
        X1, y1 = next(T_batches)


        X_adv = np.concatenate([X0, X1])
        y_class = np.concatenate([y0, np.zeros_like(y0)])

        adv_weights = []
        for layer in model.layers:
            if (layer.name.startswith("do")):
                adv_weights.append(layer.get_weights())

        if(enable_dann):
            # note - even though we save and append weights, the batchnorms moving means and variances
            # are not saved throught this mechanism 
            model.train_on_batch(X_adv, [y_class, y_adversarial_1],
                                     sample_weight=[sample_weights_class, sample_weights_adversarial])
            
            k = 0
            for layer in model.layers:
                if (layer.name.startswith("do")):
                    layer.set_weights(adv_weights[k])
                    k += 1

            class_weights = []
            
        
            for layer in model.layers:
                if (not layer.name.startswith("do")):
                    class_weights.append(layer.get_weights())
            
            domain_classification_model.train_on_batch(X_adv, [y_adversarial_2])

            k = 0
            for layer in model.layers:
                if (not layer.name.startswith("do")):
                    layer.set_weights(class_weights[k])
                    k += 1

        else:
            source_classification_model.train_on_batch(X0,y0)
            
        
        if yt is None:
            if ((i + 1) % 100 == 0):
                # print(i, stats)
                y_test_hat_s = source_classification_model.predict(Xs).argmax(1)
                domainloss,domainacc  = domain_classification_model.evaluate(np.concatenate([Xs, Xt]),
                                                                     to_categorical(np.array(([1] * Xs.shape[0] + [0] * Xt.shape[0]))),
                                                                     verbose=0)
                print("Iteration %d, source accuracy =  %.3f, discriminator loss = %.3f"%(i, accuracy_score(ys, y_test_hat_s),domainacc))
        else:
            if ((i + 1) % 100 == 0):
                # print(i, stats)
                y_test_hat_t = source_classification_model.predict(Xt).argmax(1)
                y_test_hat_s = source_classification_model.predict(Xs).argmax(1)
                print("Iteration %d, source accuracy =  %.3f, target accuracy = %.3f"%(i, accuracy_score(ys, y_test_hat_s), accuracy_score(yt, y_test_hat_t)))
                
    return embeddings_model, source_classification_model 

if __name__ == '__main__':
    Xs, ys = make_blobs(300, centers=[[0, 0], [0, 1], [0,2]], cluster_std=0.2)
    Xt, yt = make_blobs(300, centers=[[1, -1], [1.5, 0],[2,1]], cluster_std=0.2)
    
#    Xs, ys = make_blobs(300, n_features = 10, centers=None, cluster_std=0.2, 
#                        center_box = (-1,1), random_state=10)
#    Xt, yt= make_blobs(300, n_features = 10, centers=None, cluster_std=0.2, 
#                        center_box = (-1,1), random_state=3)
    
    plt.scatter(Xs[:,0], Xs[:,1], c=[["lightgreen", "tomato","lightblue"][k] for k in ys], alpha=0.4)
    plt.scatter(Xt[:,0], Xt[:,1], c=[["g", "red","blue"][k] for k in yt],  alpha=0.4)
    
    
    embs, clssmodel = train(Xs, ys, Xt, enable_dann = True, alpha=10, alpha_lr=100)
    Xs_act = embs.predict(Xs)
    Xt_act = embs.predict(Xt)
    
    plt.figure()
    plt.scatter(Xs_act[:,0], Xs_act[:,1], c=[["lightgreen", "tomato","lightblue"][k] for k in ys], alpha=0.4)
    plt.scatter(Xt_act[:,0], Xt_act[:,1], c=[["g", "red","blue"][k] for k in yt],  alpha=0.4)