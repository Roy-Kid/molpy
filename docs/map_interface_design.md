# 新的Map接口设计总结

## 🎯 设计理念

不再在 `Compute` 类里写死 `SelectResult`，而是提供一个灵活的 `map` 方法，允许用户从多种来源获取要映射的实体ID。

## 🔧 核心接口

### Compute.map() 方法

```python
def map(self, context: ComputeContext, 
        entity_ids: Optional[np.ndarray] = None,
        selection_field: str = "molecules") -> ComputeContext:
```

**参数:**
- `context`: 输入上下文
- `entity_ids`: 可选的实体ID数组。如果为None，则使用selection_field中的所有实体
- `selection_field`: 当entity_ids为None时，从哪个字段获取实体

**返回:**
- 包含映射结果的新上下文，结果存储在:
  - `result["{field}_results"]`: 字典，key为实体ID，value为计算结果
  - `result["mapped_{field}"]`: 被映射的实体ID数组

## 🛠️ 工具函数

### get_selected_entities()
```python
def get_selected_entities(context: ComputeContext) -> Optional[np.ndarray]:
```
从上下文链中提取已选择的实体ID（如果有SelectResult的话）。

### get_entities_from_frame()
```python
def get_entities_from_frame(context: ComputeContext, field: str) -> Optional[np.ndarray]:
```
从frame的指定字段获取所有实体ID。

## 📝 使用示例

### 1. 使用显式实体ID
```python
gyration_calc = GyrationCompute("gyration")

# 直接指定要计算的分子ID
specific_ids = np.array([0, 2, 4])
result = gyration_calc.map(context, entity_ids=specific_ids, selection_field="molecules")
```

### 2. 使用选择结果
```python
# 先进行选择
selection = ExpressionSelection("even_mols", select_even_func, select_field="molecules")
selected_context = selection.compute(context)

# 获取选择的ID
selected_ids = get_selected_entities(selected_context)

# 映射计算
result = gyration_calc.map(context, entity_ids=selected_ids, selection_field="molecules")
```

### 3. 映射所有实体
```python
# 计算所有分子的回转半径
result = gyration_calc.map(context, selection_field="molecules")
```

### 4. 使用工具函数
```python
# 获取所有分子ID
all_mol_ids = get_entities_from_frame(context, "molecules")

# 选择子集（比如每隔一个分子）
subset_ids = all_mol_ids[::2]

# 映射计算
result = gyration_calc.map(context, entity_ids=subset_ids, selection_field="molecules")
```

## ✅ 优势

1. **灵活性**: 不依赖特定的选择机制，entity_ids可以来自任何地方
2. **解耦**: Compute类不需要知道SelectResult的存在
3. **可组合**: 可以轻松组合不同的选择逻辑
4. **工具丰富**: 提供工具函数帮助从不同来源获取entity_ids
5. **类型安全**: 清晰的接口和类型注解

## 🧪 测试覆盖

- ✅ 显式entity_ids映射
- ✅ 全实体映射
- ✅ 与SelectResult结合使用
- ✅ 工具函数功能
- ✅ 不同字段映射（分子、原子等）
- ✅ 错误情况处理

## 📊 使用场景

1. **分子性质计算**: 对选定分子计算回转半径、偶极矩等
2. **原子分析**: 对特定原子进行电荷、坐标分析
3. **团簇研究**: 对分子团簇进行各种几何或能量计算
4. **时间序列**: 对轨迹中的特定帧进行分析

这个设计完全符合你的要求：
- ❌ 不在Compute类里写死SelectResult
- ❌ 没有VectorizedResult和VectorizedCompute
- ✅ 提供了简洁的map方法实现map-reduce模式
- ✅ 完全灵活，支持任意来源的entity_ids
