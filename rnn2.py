import tensorflow as tf
from tensorflow.models.rnn import rnn, rnn_cell
from sklearn.model_selection import KFold
import numpy as np
from pyfasta import Fasta
from util import *
import sys
import random

seqs = Fasta('data/hg19.fa')
pos = read_bed('data/hsa_hg19_Rybak2015.bed')
neg = read_bed('tmp/negatives.bed')
data = pos + neg
labels = np.load('tmp/labels.npy')
test, train = list(KFold(n_splits=2, shuffle=True).split(data))[0]
print("Data loaded")

# Parameters
learning_rate = 0.001
training_iters = 100000
batch_size = 128
display_step = 1

# Network Parameters
n_input = 5 # Input (ACGTN one hot)
n_steps = 100 #max(g.end - g.start for g in data) # timesteps
n_hidden = 64 # hidden layer num of features
n_classes = 2

# tf Graph input
x = tf.placeholder("float", [None, n_steps, n_input])
istate = tf.placeholder("float", [None, n_hidden]) #state & cell => 2x n_hidden
y = tf.placeholder("float", [None, n_classes])

# Define weights
weights = {
	'hidden': tf.Variable(tf.random_normal([n_input, n_hidden])), # Hidden layer weights
	'out': tf.Variable(tf.random_normal([n_hidden, n_classes]))
}
biases = {
	'hidden': tf.Variable(tf.random_normal([n_hidden])),
	'out': tf.Variable(tf.random_normal([n_classes]))
}

def RNN(_X, _istate, _weights, _biases):
	# input shape: (batch_size, n_steps, n_input)
	_X = tf.transpose(_X, [1, 0, 2])  # permute n_steps and batch_size
	# Reshape to prepare input to hidden activation
	_X = tf.reshape(_X, [-1, n_input]) # (n_steps*batch_size, n_input)
	# Linear activation
	_X = tf.matmul(_X, _weights['hidden']) + _biases['hidden']

	# Define a lstm cell with tensorflow
	lstm_cell = rnn_cell.GRUCell(n_hidden)
	# Split data because rnn cell needs a list of inputs for the RNN inner loop
	_X = tf.split(0, n_steps, _X) # n_steps * (batch_size, n_hidden)

	# Get lstm cell output
	outputs, states = rnn.rnn(lstm_cell, _X, initial_state=_istate)

	# Linear activation
	# Get inner loop last output
	return tf.matmul(outputs[-1], _weights['out']) + _biases['out']

pred = RNN(x, istate, weights, biases)

# Define loss and optimizer
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(pred, y)) # Softmax loss
optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost) # Adam Optimizer

# Evaluate model
correct_pred = tf.equal(tf.argmax(pred,1), tf.argmax(y,1))
accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32))

# Initializing the variables
init = tf.initialize_all_variables()

print("Graph built, steps: %d" % n_steps)
sys.stdout.flush()

# Launch the graph
with tf.Session() as sess:
	sess.run(init)
	step = 1
	# Keep training until reach max iterations
	while step * batch_size < training_iters:
		batch = random.sample(list(train), batch_size)
		batch_xs, batch_ys = data[batch], labels[batch]

		# Get actual DNA data
		batch_xs = [ 
			get_gene_data(g, seqs, one_hot=True, pad_len=n_steps)[:n_steps] 
			for g in batch_xs ]

		# One-hot encoding for classifier
		batch_ys = [ [ x, 1-x ] for x in batch_ys ]

		# Reshape data to get 28 seq of 28 elements
		batch_xs = batch_xs.reshape((batch_size, n_steps, n_input))
		# Fit training using batch data
		sess.run(optimizer, feed_dict={x: batch_xs, y: batch_ys,
									   istate: np.zeros((batch_size, 2*n_hidden))})
		if step % display_step == 0:
			# Calculate batch accuracy
			acc = sess.run(accuracy, feed_dict={x: batch_xs, y: batch_ys,
												istate: np.zeros((batch_size, 2*n_hidden))})
			# Calculate batch loss
			loss = sess.run(cost, feed_dict={x: batch_xs, y: batch_ys,
											 istate: np.zeros((batch_size, 2*n_hidden))})
			print("Iter " + str(step*batch_size) + ", Minibatch Loss= " + "{:.6f}".format(loss) + \
				  ", Training Accuracy= " + "{:.5f}".format(acc))
			sys.stdout.flush()
		step += 1
	print("Optimization Finished!")
	# Calculate accuracy for 256 mnist test images
	test_len = 256
	test_data = mnist.test.images[:test_len].reshape((-1, n_steps, n_input))
	test_label = mnist.test.labels[:test_len]
	print("Testing Accuracy:", sess.run(accuracy, feed_dict={x: test_data, y: test_label,
															 istate: np.zeros((test_len, 2*n_hidden))}))