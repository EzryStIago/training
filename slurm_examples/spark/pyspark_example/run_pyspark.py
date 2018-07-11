
from pyspark.sql import SparkSession
from pyspark.sql.functions import *

# One way of chaining methods in Python is () around the multiline expression
spark_builder = (
    SparkSession
    .builder
    .appName("GroupIrises"))

# SparkSession is the entry point in the application
spark_sess = spark_builder.getOrCreate()

# read a csv file which has a header and infer the schema from the data (very slow!)
# Another way to chain methods in python is line continuation character \ with nothing after it
df = spark_sess.read                                             \
               .option("header", "true")                         \
               .option("inferSchema", "true")                    \
               .csv("input_data/iris.csv")

# group by column species, compute the average per group of sepal_length, 
# the resulting dataframe might have lots of partitions (chunks) 
# merge them all to 1 since it is a small dataset, with coalesce
# write out the result with headers and overwrite any previous file in that spot
df.groupBy("species")                                           \
    .mean("sepal_length")                                       \
    .coalesce(1)                                                \
    .write.option("header", "true")                             \
    .option("mode", "overwrite")                                \
    .csv("output_data/grouped_iris.csv")

spark_sess.stop()


