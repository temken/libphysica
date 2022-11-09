#ifndef __List_Manipulations_hpp_
#define __List_Manipulations_hpp_

namespace libphysica
{

template <typename T>
extern bool Lists_Equal(const std::vector<T>& v1, const std::vector<T>& v2)
{
	if(v1.size() != v2.size())
		return false;
	for(unsigned int i = 0; i < v1.size(); i++)
		if(v1[i] != v2[i])
			return false;
	return true;
}

template <typename T>
extern bool Lists_Equal(const std::vector<std::vector<T>>& v1, const std::vector<std::vector<T>>& v2)
{
	if(v1.size() != v2.size())
		return false;
	for(unsigned int i = 0; i < v1.size(); i++)
		if(!Lists_Equal(v1[i], v2[i]))
			return false;
	return true;
}

template <typename T>
extern std::vector<T> Combine_Lists(const std::vector<T>& v1, const std::vector<T>& v2)
{
	std::vector<T> result = v1;
	result.insert(result.end(), v2.begin(), v2.end());
	return result;
}

template <typename T>
extern std::vector<std::vector<T>> Transpose_Lists(const std::vector<T>& v1, const std::vector<T>& v2)
{
	if(v1.size() != v2.size())
	{
		std::cerr << "Error in libphysica::Transpose_Vectors(): Vectors have different sizes." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::vector<std::vector<T>> result;
	for(unsigned int i = 0; i < v1.size(); i++)
		result.push_back({v1[i], v2[i]});
	return result;
}

template <typename T>
extern std::vector<T> Sub_List(const std::vector<T>& v, int i1, int i2)
{
	if(i1 < 0)
		i1 = 0;
	if(i2 > v.size())
		i2 = v.size();
	std::vector<T> sub(&v[i1], &v[i2] + 1);
	return sub;
}

template <typename T>
extern std::vector<T> Flatten_List(const std::vector<std::vector<T>>& v)
{
	std::vector<T> result;
	for(unsigned int i = 0; i < v.size(); i++)
		result.insert(result.end(), v[i].begin(), v[i].end());
	return result;
}

template <typename T>
extern bool List_Contains(const std::vector<T>& list, T value)
{
	for(unsigned int i = 0; i < list.size(); i++)
		if(list[i] == value)
			return true;
	return false;
}

template <typename T>
extern std::vector<int> Find_Indices(const std::vector<T>& v, T value)
{
	std::vector<int> indices;
	for(unsigned int i = 0; i < v.size(); i++)
		if(v[i] == value)
			indices.push_back(i);
	return indices;
}

}	// namespace libphysica

#endif